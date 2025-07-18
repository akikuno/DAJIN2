from __future__ import annotations

import csv
import json
import logging
import os
import queue
import shutil
import socket
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


def upload_files(files, expected_type=None):
    """Upload files with validation"""
    for file in files:
        # Validate file
        file_type = validate_file(file)
        if expected_type and file_type not in expected_type:
            raise ValueError(f"Unexpected file format. Expected {expected_type} but got {file_type}")

        path = Path(app.config["UPLOAD_FOLDER"], secure_filename(file.filename))
        file.save(path)


def get_path_of_uploaded_files(UPLOAD_FOLDER, files):
    paths = []
    for file in files:
        path = str(UPLOAD_FOLDER) + "/" + secure_filename(file.filename)
        paths.append(path)
    return paths


app = Flask(__name__)


@app.route("/")
def root_page():
    return render_template("index.html")


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
        threads = request.form.get("batch-threads")
        try:
            threads = int(threads) if threads else 1
            if threads < 1:
                threads = 1
        except ValueError:
            threads = 1

        no_filter = request.form.get("batch-no-filter") == "on"

        # Prepare arguments
        arguments = {"file": str(batch_file_path), "threads": threads, "debug": False, "no_filter": no_filter}

        # Create progress queue for this analysis
        progress_queue = queue.Queue()
        progress_queues[batch_id] = progress_queue

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
                progress_queue.put(
                    {"status": "log", "message": "Starting batch analysis...", "timestamp": time.time()}
                )

                # Run the actual batch analysis
                main.execute_batch_mode(arguments)

                # Store result
                analysis_results[batch_id] = {
                    "status": "completed",
                    "result_path": "DAJIN_Results",
                    "abs_result_path": str(Path("DAJIN_Results").resolve()),
                    "timestamp": time.time(),
                }

                completion_message = (
                    f"🎊 Batch analysis completed! Your results are saved in {str(Path('DAJIN_Results').resolve())}"
                )
                progress_queue.put(
                    {
                        "status": "completed",
                        "message": completion_message,
                        "result_path": str(Path("DAJIN_Results").resolve()),
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
                except:
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
                except:
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
        name = request.form.get("name")
        if not name or not name.strip():
            return jsonify({"error": "Please enter a project name"}), 400

        name = name.strip()

        # Validate required files
        if "sample" not in request.files or not request.files.getlist("sample"):
            return jsonify({"error": "Please select sample files"}), 400

        if "control" not in request.files or not request.files.get("control").filename:
            return jsonify({"error": "Please select a control file"}), 400

        if "allele" not in request.files or not request.files.getlist("allele"):
            return jsonify({"error": "Please select allele files"}), 400

        # Setup directories
        TEMPDIR = Path(config.TEMP_ROOT_DIR, name)
        if TEMPDIR.exists():
            shutil.rmtree(TEMPDIR)
        UPLOAD_FOLDER = Path(TEMPDIR, "upload")
        UPLOAD_FOLDER.mkdir(parents=True, exist_ok=True)
        app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER

        # Handle directory structure for DAJIN2
        # DAJIN2 expects directories for sample/control, files for allele

        # Process sample directory upload
        sample_files = request.files.getlist("sample")
        if not sample_files:
            return jsonify({"error": "Please select sample directory"}), 400

        sample_dir = UPLOAD_FOLDER / "sample"
        sample_dir.mkdir(exist_ok=True)

        # Preserve directory structure from webkitdirectory upload
        for file in sample_files:
            if file.filename:
                # Validate file type
                validate_file(file)
                # Get relative path from the uploaded file
                relative_path = file.filename
                # Create directory structure
                file_path = sample_dir / Path(relative_path).name  # Just use filename, ignore subdirs for now
                file.save(file_path)

        PATH_SAMPLE = str(sample_dir)

        # Process control directory upload
        control_files = request.files.getlist("control")
        if not control_files:
            return jsonify({"error": "Please select control directory"}), 400

        control_dir = UPLOAD_FOLDER / "control"
        control_dir.mkdir(exist_ok=True)

        for file in control_files:
            if file.filename:
                # Validate file type
                validate_file(file)
                relative_path = file.filename
                file_path = control_dir / Path(relative_path).name
                file.save(file_path)

        PATH_CONTROL = str(control_dir)

        # Upload allele files directly (DAJIN2 expects file path for alleles)
        allele_files = request.files.getlist("allele")
        if not allele_files:
            return jsonify({"error": "Please select allele files"}), 400

        upload_files(allele_files, expected_type=["fasta"])
        PATH_ALLELE = get_path_of_uploaded_files(UPLOAD_FOLDER, allele_files)

        # Handle BED file upload (optional)
        PATH_BED = None
        if "bed" in request.files and request.files.get("bed").filename:
            bed_file = request.files["bed"]
            file_type = validate_file(bed_file)
            if file_type != "bed":
                return jsonify({"error": "BED file must have .bed extension"}), 400
            bed_path = Path(UPLOAD_FOLDER, secure_filename(bed_file.filename))
            bed_file.save(bed_path)
            PATH_BED = str(bed_path)

        genome = request.form.get("genome", "").strip()
        threads = request.form.get("threads")
        try:
            threads = int(threads) if threads else 1
            if threads < 1:
                threads = 1
        except ValueError:
            threads = 1

        # Get no-filter option
        no_filter = request.form.get("no-filter") == "on"

        # Prepare batch data
        data = []
        row = {"sample": PATH_SAMPLE, "name": name, "control": PATH_CONTROL, "allele": PATH_ALLELE[0]}
        if genome:
            row["genome"] = genome
        if PATH_BED:
            row["bed"] = PATH_BED
        data.append(row)

        # Create batch file
        batch_file_path = Path(TEMPDIR, "upload", "batch.csv")
        batch_file_path.parent.mkdir(parents=True, exist_ok=True)

        with open(batch_file_path, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(data)

        # Prepare arguments
        arguments = {"file": str(batch_file_path), "threads": threads, "debug": False, "no_filter": no_filter}

        # Create progress queue for this analysis
        progress_queue = queue.Queue()
        progress_queues[name] = progress_queue

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

                completion_message = f"🎊 Analysis completed! Your results are saved in {abs_result_path}"
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
                except:
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
                except:
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
                        completion_message = f"🎊 Analysis completed! Your results are saved in {abs_path}"
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
    webbrowser.open_new(f"http://localhost:{PORT}/")


def execute():
    PORT = find_free_port()
    print(f"Assess 'http://localhost:{PORT}/' if a browser does not automatically open.")
    Timer(1, open_browser, [PORT]).start()
    with redirect_stderr(open(os.devnull, "w")):
        serve(app, host="0.0.0.0", port=PORT)
