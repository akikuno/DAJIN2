from __future__ import annotations

import csv
import json
import logging
import os
import queue
import shutil
import socket
import subprocess
import threading
import time
import webbrowser
from contextlib import closing, redirect_stderr
from pathlib import Path
from threading import Timer

from flask import Flask, Response, jsonify, render_template, request
from waitress import serve
from werkzeug.utils import secure_filename

from DAJIN2 import main
from DAJIN2.utils import config

# Global variable to store progress information
progress_queues = {}
analysis_results = {}
GUI_UPLOAD_ROOT_DIR = Path(config.DAJIN_RESULTS_DIR, ".gui_upload")
MIN_THREADS = 1
MAX_THREADS = 32
SINGLE_MODE_INPUT_TYPES = {"fastq", "fasta", "bam"}


def build_completion_message(is_batch: bool = False) -> str:
    scope = "batch analysis" if is_batch else "analysis"
    return f"Your {scope} results are saved in the following directory:"


def normalize_project_name(project_name: str) -> str:
    sanitized_name = secure_filename(project_name.strip())
    if not sanitized_name:
        raise ValueError("Project name must include letters, numbers, underscores, hyphens, or dots.")
    return sanitized_name


def parse_threads(raw_value: str | None, default: int = 1) -> int:
    try:
        threads = int(raw_value) if raw_value else default
    except (TypeError, ValueError):
        return default

    return max(MIN_THREADS, min(MAX_THREADS, threads))


class ProgressLogHandler(logging.Handler):
    """Custom logging handler to capture DAJIN2 log messages and send them to progress queue."""

    def __init__(self, progress_queue):
        super().__init__()
        self.progress_queue = progress_queue
        self.setLevel(logging.INFO)
        # Format to match the log file format: timestamp, level, message
        formatter = logging.Formatter("%(asctime)s, %(levelname)s, %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        self.setFormatter(formatter)

    def emit(self, record):
        try:
            # Only capture DAJIN2 related logs to avoid spam
            if record.name.startswith("DAJIN2") or "DAJIN" in record.getMessage():
                # Format the log message
                log_message = self.format(record)
                # Send to progress queue
                self.progress_queue.put({"status": "log", "message": log_message, "timestamp": time.time()})
        except Exception:
            # Don't let logging errors break the analysis
            pass


def convert_to_windows_path(unix_path):
    """Convert Unix/WSL path to Windows path"""
    import re

    # Convert WSL mount paths (/mnt/drive/) to Windows drive letters
    wsl_pattern = r"^/mnt/([a-z])(/.*)?$"
    match = re.match(wsl_pattern, unix_path.lower())

    if match:
        drive_letter = match.group(1).upper()
        rest_of_path = match.group(2) or ""
        # Convert forward slashes to backslashes
        windows_path = f"{drive_letter}:{rest_of_path}".replace("/", "\\")
        return windows_path
    else:
        # For non-WSL paths, just convert slashes
        return unix_path.replace("/", "\\")


def find_free_port():
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


def validate_file(file):
    """Validate uploaded file"""
    if not file.filename:
        raise ValueError("No file selected")

    # Check file extension
    allowed_extensions = {
        "fastq": [".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        "fasta": [".fasta", ".fa", ".fasta.gz", ".fa.gz"],
        "bam": [".bam"],
        "csv": [".csv"],
        "excel": [".xlsx", ".xls"],
        "bed": [".bed"],
    }

    filename = file.filename.lower()
    for file_type, extensions in allowed_extensions.items():
        if any(filename.endswith(ext) for ext in extensions):
            return file_type

    raise ValueError(f"Unsupported file format: {file.filename}")


def collect_valid_uploaded_files(files, expected_types: set[str], field_name: str):
    """Collect uploaded files and validate their formats."""
    valid_files = []
    invalid_files = []

    for file in files:
        if not file or not file.filename:
            continue

        try:
            file_type = validate_file(file)
        except ValueError:
            invalid_files.append(Path(file.filename).name)
            continue

        if file_type not in expected_types:
            invalid_files.append(Path(file.filename).name)
            continue

        valid_files.append(file)

    if invalid_files:
        invalid_file_names = ", ".join(invalid_files[:5])
        if len(invalid_files) > 5:
            invalid_file_names += ", ..."
        raise ValueError(f"{field_name} contains unsupported files: {invalid_file_names}")

    if not valid_files:
        raise ValueError(f"{field_name} does not contain supported files.")

    return valid_files


def save_uploaded_files_to_directory(files, output_dir: Path) -> list[str]:
    """Save uploaded files as flat files to match DAJIN2 directory expectations."""
    saved_paths = []
    for file in files:
        safe_filename = secure_filename(file.filename)
        if not safe_filename:
            continue

        output_path = output_dir / safe_filename
        file.save(output_path)
        saved_paths.append(str(output_path))

    return saved_paths


app = Flask(__name__)


def is_wsl_environment() -> bool:
    if os.name != "posix":
        return False

    if os.environ.get("WSL_DISTRO_NAME") or os.environ.get("WSL_INTEROP"):
        return True

    try:
        with open("/proc/sys/kernel/osrelease") as kernel_info:
            return "microsoft" in kernel_info.read().lower()
    except OSError:
        return False


def open_url_in_windows_default_browser(url: str) -> bool:
    commands = [
        ["cmd.exe", "/c", "start", "", url],
        ["powershell.exe", "-NoProfile", "-Command", f'Start-Process "{url}"'],
    ]

    for command in commands:
        try:
            subprocess.run(
                command,
                check=True,
                timeout=5,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return True
        except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
            continue

    return False


@app.route("/")
def root_page():
    return render_template("gui.html")


@app.route("/submit-batch", methods=["POST"])
def submit_batch():
    try:
        # Validate batch file
        if "batch-file" not in request.files or not request.files.get("batch-file").filename:
            return jsonify({"error": "Please select a batch file"}), 400

        batch_file = request.files["batch-file"]
        file_type = validate_file(batch_file)

        if file_type not in ["csv", "excel"]:
            return jsonify({"error": "Batch file must be CSV or Excel format"}), 400

        # Setup directories
        import uuid

        batch_id = f"batch_{uuid.uuid4().hex[:8]}"
        TEMPDIR = Path(config.TEMP_ROOT_DIR, batch_id)
        if TEMPDIR.exists():
            shutil.rmtree(TEMPDIR)
        UPLOAD_FOLDER = Path(TEMPDIR, "upload")
        UPLOAD_FOLDER.mkdir(parents=True, exist_ok=True)

        # Save batch file
        batch_filename = secure_filename(batch_file.filename)
        batch_file_path = UPLOAD_FOLDER / batch_filename
        batch_file.save(batch_file_path)

        # Get threads and no-filter options
        threads = parse_threads(request.form.get("batch-threads"))

        no_filter = request.form.get("batch-no-filter") == "on"

        # Create progress queue for this analysis
        progress_queue = queue.Queue()
        progress_queues[batch_id] = progress_queue

        # Prepare arguments
        arguments = {
            "file": str(batch_file_path),
            "threads": threads,
            "debug": False,
            "no_filter": no_filter,
            "progress_queue": progress_queue,
        }

        # Start analysis in background thread
        def run_batch_analysis():
            try:
                # Set up log handler to capture DAJIN2 log messages
                log_handler = ProgressLogHandler(progress_queue)

                # Add handler to root logger to capture ALL logs (filtered in handler)
                root_logger = logging.getLogger()
                root_logger.addHandler(log_handler)

                # Also add to DAJIN2 main logger for good measure
                dajin_logger = logging.getLogger("DAJIN2")
                dajin_logger.addHandler(log_handler)
                dajin_logger.setLevel(logging.INFO)

                # Send initial status
                progress_queue.put({"status": "log", "message": "Starting batch analysis...", "timestamp": time.time()})

                # Run the actual batch analysis
                main.execute_batch_mode(arguments)

                # Store result
                analysis_results[batch_id] = {
                    "status": "completed",
                    "result_path": "DAJIN_Results",
                    "abs_result_path": str(Path("DAJIN_Results").resolve()),
                    "timestamp": time.time(),
                }

                result_directory = str(Path("DAJIN_Results").resolve())
                completion_message = build_completion_message(is_batch=True)
                progress_queue.put(
                    {
                        "status": "completed",
                        "message": completion_message,
                        "result_path": result_directory,
                    }
                )

                # Clean up handlers
                root_logger.removeHandler(log_handler)
                dajin_logger.removeHandler(log_handler)

            except FileNotFoundError as e:
                # Clean up handlers before handling error
                try:
                    root_logger.removeHandler(log_handler)
                    dajin_logger.removeHandler(log_handler)
                except (UnboundLocalError, AttributeError):
                    pass

                # Enhanced error message for file/directory not found
                current_dir = Path.cwd().resolve()
                error_str = str(e)

                # Extract the missing path from the error message
                import re

                path_match = re.search(r"'([^']*)'", error_str)
                missing_path = path_match.group(1) if path_match else "unknown path"

                detailed_error = (
                    f"File or directory not found: '{missing_path}'\n\n"
                    f"Current working directory: {current_dir}\n\n"
                    f"The path '{missing_path}' does not exist in the current directory. "
                    f"Please check the following:\n"
                    f"1. Update the paths in your batch file to match the actual file locations\n"
                    f"2. Ensure all sample, control, and allele files are accessible from '{current_dir}'\n"
                    f"3. Use absolute paths in your batch file, or\n"
                    f"4. Change the working directory where you started DAJIN2 GUI to match your batch file paths\n\n"
                    f"Original error: {error_str}"
                )

                progress_queue.put({"status": "error", "message": detailed_error})
                analysis_results[batch_id] = {"status": "error", "error": detailed_error, "timestamp": time.time()}

            except Exception as e:
                # Clean up handlers before handling error
                try:
                    root_logger.removeHandler(log_handler)
                    dajin_logger.removeHandler(log_handler)
                except (UnboundLocalError, AttributeError):
                    pass

                error_msg = f"Batch analysis error occurred: {str(e)}"
                progress_queue.put({"status": "error", "message": error_msg})
                analysis_results[batch_id] = {"status": "error", "error": error_msg, "timestamp": time.time()}

        # Start analysis thread
        analysis_thread = threading.Thread(target=run_batch_analysis)
        analysis_thread.daemon = True
        analysis_thread.start()

        return jsonify({"status": "started", "message": "Batch analysis started", "analysis_id": batch_id})

    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": f"Unexpected error occurred: {str(e)}"}), 500


@app.route("/submit", methods=["POST"])
def submit():
    try:
        # Get form data
        name = request.form.get("name", "")
        if not name.strip():
            return jsonify({"error": "Please enter a project name"}), 400

        name = normalize_project_name(name)

        # Validate required files
        sample_files = request.files.getlist("sample")
        if not sample_files:
            return jsonify({"error": "Please select sample files"}), 400

        control_files = request.files.getlist("control")
        if not control_files:
            return jsonify({"error": "Please select a control file"}), 400

        allele_files = request.files.getlist("allele")
        if not allele_files:
            return jsonify({"error": "Please select allele files"}), 400

        sample_files = collect_valid_uploaded_files(sample_files, SINGLE_MODE_INPUT_TYPES, "Sample directory")
        control_files = collect_valid_uploaded_files(control_files, SINGLE_MODE_INPUT_TYPES, "Control directory")
        allele_files = collect_valid_uploaded_files(allele_files, {"fasta"}, "Allele file")

        if len(allele_files) != 1:
            return jsonify({"error": "Please select exactly one allele FASTA file (.fa or .fasta)"}), 400

        # Setup directories.
        # Keep GUI uploads outside DAJIN2 runtime tempdir to avoid cleanup collisions.
        GUI_TEMPDIR = Path(GUI_UPLOAD_ROOT_DIR, name)
        if GUI_TEMPDIR.exists():
            shutil.rmtree(GUI_TEMPDIR)
        UPLOAD_FOLDER = Path(GUI_TEMPDIR, "upload")
        UPLOAD_FOLDER.mkdir(parents=True, exist_ok=True)
        app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER

        # Handle directory structure for DAJIN2
        # DAJIN2 expects directories for sample/control, files for allele

        # Process sample directory upload
        sample_dir = UPLOAD_FOLDER / "sample"
        sample_dir.mkdir(parents=True, exist_ok=True)
        save_uploaded_files_to_directory(sample_files, sample_dir)

        PATH_SAMPLE = str(sample_dir)

        # Process control directory upload
        control_dir = UPLOAD_FOLDER / "control"
        control_dir.mkdir(parents=True, exist_ok=True)
        save_uploaded_files_to_directory(control_files, control_dir)

        PATH_CONTROL = str(control_dir)

        # Upload a single allele file directly (DAJIN2 expects one path in the batch CSV).
        allele_file = allele_files[0]
        allele_name = secure_filename(Path(allele_file.filename).name)
        allele_path = UPLOAD_FOLDER / allele_name
        allele_file.save(allele_path)
        PATH_ALLELE = str(allele_path)

        # Handle BED file upload (optional)
        PATH_BED = None
        if "bed" in request.files and request.files.get("bed").filename:
            bed_file = request.files["bed"]
            bed_file = collect_valid_uploaded_files([bed_file], {"bed"}, "BED file")[0]
            bed_name = secure_filename(Path(bed_file.filename).name)
            bed_path = Path(UPLOAD_FOLDER, bed_name)
            bed_file.save(bed_path)
            PATH_BED = str(bed_path)

        genome = request.form.get("genome", "").strip()
        threads = parse_threads(request.form.get("threads"))

        # Get no-filter option
        no_filter = request.form.get("no-filter") == "on"

        # Prepare batch data
        data = []
        row = {"sample": PATH_SAMPLE, "name": name, "control": PATH_CONTROL, "allele": PATH_ALLELE}
        if genome:
            row["genome"] = genome
        if PATH_BED:
            row["bed"] = PATH_BED
        data.append(row)

        # Create batch file
        batch_file_path = Path(GUI_TEMPDIR, "upload", "batch.csv")
        batch_file_path.parent.mkdir(parents=True, exist_ok=True)

        with open(batch_file_path, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(data)

        # Create progress queue for this analysis
        progress_queue = queue.Queue()
        progress_queues[name] = progress_queue

        # Prepare arguments
        arguments = {
            "file": str(batch_file_path),
            "threads": threads,
            "debug": False,
            "no_filter": no_filter,
            "progress_queue": progress_queue,
        }

        # Start analysis in background thread
        def run_analysis():
            try:
                # Set up log handler to capture DAJIN2 log messages
                log_handler = ProgressLogHandler(progress_queue)

                # Add handler to root logger to capture ALL logs (filtered in handler)
                root_logger = logging.getLogger()
                root_logger.addHandler(log_handler)

                # Also add to DAJIN2 main logger for good measure
                dajin_logger = logging.getLogger("DAJIN2")
                dajin_logger.addHandler(log_handler)
                dajin_logger.setLevel(logging.INFO)

                # Send initial status
                progress_queue.put({"status": "log", "message": "Starting analysis...", "timestamp": time.time()})

                # Run the actual analysis
                main.execute_batch_mode(arguments)

                # Store result
                result_path = f"DAJIN_Results/{name}"
                abs_result_path = str(Path(result_path).resolve())
                analysis_results[name] = {
                    "status": "completed",
                    "result_path": result_path,
                    "abs_result_path": abs_result_path,
                    "timestamp": time.time(),
                }

                completion_message = build_completion_message(abs_result_path)
                progress_queue.put(
                    {"status": "completed", "message": completion_message, "result_path": abs_result_path}
                )

                # Clean up handlers
                root_logger.removeHandler(log_handler)
                dajin_logger.removeHandler(log_handler)

            except FileNotFoundError as e:
                # Clean up handlers before handling error
                try:
                    root_logger.removeHandler(log_handler)
                    dajin_logger.removeHandler(log_handler)
                except (UnboundLocalError, AttributeError):
                    pass

                # Enhanced error message for file/directory not found
                current_dir = Path.cwd().resolve()
                error_str = str(e)

                # Extract the missing path from the error message
                import re

                path_match = re.search(r"'([^']*)'", error_str)
                missing_path = path_match.group(1) if path_match else "unknown path"

                detailed_error = (
                    f"File or directory not found: '{missing_path}'\n\n"
                    f"Current working directory: {current_dir}\n\n"
                    f"The path '{missing_path}' does not exist in the current directory. "
                    f"Please check the following:\n"
                    f"1. Ensure all uploaded files were saved correctly\n"
                    f"2. Verify that sample and control directories contain the expected files\n"
                    f"3. Check that allele and BED files (if provided) are accessible\n"
                    f"4. Try uploading the files again if the issue persists\n\n"
                    f"Original error: {error_str}"
                )

                progress_queue.put({"status": "error", "message": detailed_error})
                analysis_results[name] = {"status": "error", "error": detailed_error, "timestamp": time.time()}

            except Exception as e:
                # Clean up handlers before handling error
                try:
                    root_logger.removeHandler(log_handler)
                    dajin_logger.removeHandler(log_handler)
                except (UnboundLocalError, AttributeError):
                    pass

                error_msg = f"Analysis error occurred: {str(e)}"
                progress_queue.put({"status": "error", "message": error_msg})
                analysis_results[name] = {"status": "error", "error": error_msg, "timestamp": time.time()}

        # Start analysis thread
        analysis_thread = threading.Thread(target=run_analysis)
        analysis_thread.daemon = True
        analysis_thread.start()

        return jsonify({"status": "started", "message": "Analysis started", "analysis_id": name})

    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": f"Unexpected error occurred: {str(e)}"}), 500


@app.route("/progress/<analysis_id>")
def get_progress(analysis_id):
    """Server-Sent Events endpoint for progress updates"""

    def generate():
        if analysis_id not in progress_queues:
            yield f"data: {json.dumps({'status': 'error', 'message': 'Analysis not found'})}\n\n"
            return

        progress_queue = progress_queues[analysis_id]

        while True:
            try:
                # Wait for progress update with timeout
                update = progress_queue.get(timeout=1)
                yield f"data: {json.dumps(update)}\n\n"

                # If analysis is completed or failed, stop sending updates
                if update.get("status") in ["completed", "error"]:
                    break

            except queue.Empty:
                # Send heartbeat to keep connection alive
                yield f"data: {json.dumps({'status': 'heartbeat'})}\n\n"

                # Check if analysis result is available (in case we missed the queue message)
                if analysis_id in analysis_results:
                    result = analysis_results[analysis_id]
                    if result["status"] == "completed":
                        abs_path = result.get("abs_result_path", result["result_path"])
                        is_batch = analysis_id.startswith("batch_")
                        completion_message = build_completion_message(abs_path, is_batch=is_batch)
                        yield f"data: {json.dumps({'status': 'completed', 'message': completion_message, 'result_path': abs_path})}\n\n"
                        break
                    elif result["status"] == "error":
                        yield f"data: {json.dumps({'status': 'error', 'message': result['error']})}\n\n"
                        break
            except Exception as e:
                yield f"data: {json.dumps({'status': 'error', 'message': f'Progress monitoring error: {str(e)}'})}\n\n"
                break

    return Response(generate(), mimetype="text/event-stream")


@app.route("/status/<analysis_id>")
def get_status(analysis_id):
    """Get current status of analysis"""
    if analysis_id in analysis_results:
        return jsonify(analysis_results[analysis_id])
    elif analysis_id in progress_queues:
        return jsonify({"status": "running", "message": "Analysis running..."})
    else:
        return jsonify({"status": "not_found", "message": "Analysis not found"}), 404


@app.route("/open-folder/<analysis_id>")
def open_result_folder(analysis_id):
    """Open result folder in system file manager"""
    import os
    import platform
    import subprocess

    try:
        if analysis_id not in analysis_results:
            return jsonify({"error": "Analysis not found"}), 404

        result = analysis_results[analysis_id]
        if result["status"] != "completed":
            return jsonify({"error": "Analysis not completed"}), 400

        result_path = result["result_path"]
        abs_path = str(Path(result_path).resolve())

        # Check if directory exists
        if not Path(abs_path).exists():
            return jsonify({"error": f"Result directory not found: {abs_path}"}), 404

        # Open folder in system file manager
        system = platform.system()

        # Convert path for Windows display
        windows_path = convert_to_windows_path(abs_path) if system == "Windows" else abs_path

        try:
            if system == "Windows":
                # Try multiple methods for Windows
                try:
                    # Method 1: Use explorer with Windows path
                    subprocess.run(["explorer", windows_path], check=True, timeout=5)
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
                    try:
                        # Method 2: Use os.startfile with Windows path
                        os.startfile(windows_path)
                    except OSError:
                        # Method 3: Try with Unix path as fallback
                        os.startfile(abs_path)

            elif system == "Darwin":  # macOS
                subprocess.run(["open", abs_path], check=True, timeout=5)
            elif system == "Linux":
                subprocess.run(["xdg-open", abs_path], check=True, timeout=5)
            else:
                return jsonify(
                    {
                        "error": f"Unsupported operating system: {system}",
                        "path": abs_path,
                        "windows_path": windows_path,
                    }
                ), 400

        except Exception as open_error:
            # Log the specific error for debugging
            error_details = f"System: {system}, Path: {abs_path}, Error: {str(open_error)}"
            return jsonify(
                {
                    "error": f"Could not open folder: {str(open_error)}",
                    "debug_info": error_details,
                    "path": abs_path,
                    "windows_path": windows_path,
                    "suggestion": "Please manually navigate to the path shown above",
                }
            ), 500

        return jsonify(
            {
                "success": True,
                "message": "Folder opened successfully",
                "path": abs_path,
                "windows_path": windows_path,
                "system": system,
            }
        )

    except Exception as e:
        return jsonify({"error": f"Failed to open folder: {str(e)}"}), 500


def open_browser(PORT):
    url = f"http://localhost:{PORT}/"

    if is_wsl_environment():
        if not open_url_in_windows_default_browser(url):
            print(f"Could not open the host browser automatically. Access '{url}' manually.")
        return

    webbrowser.open_new(url)


def execute():
    PORT = find_free_port()
    print(f"Assess 'http://localhost:{PORT}/' if a browser does not automatically open.")
    Timer(1, open_browser, [PORT]).start()
    with open(os.devnull, "w") as devnull, redirect_stderr(devnull):
        serve(app, host="0.0.0.0", port=PORT)
