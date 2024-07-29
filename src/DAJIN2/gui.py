from __future__ import annotations

import csv
import os
import shutil
import socket
import webbrowser
from contextlib import closing, redirect_stderr
from pathlib import Path
from threading import Timer

from flask import Flask, render_template, request
from waitress import serve
from werkzeug.utils import secure_filename

from DAJIN2 import main
from DAJIN2.utils import config


def find_free_port():
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


def upload_files(files):
    for file in files:
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


@app.route("/submit", methods=["POST"])
def submit():
    name = request.form.get("name")
    TEMPDIR = Path(config.TEMP_ROOT_DIR, name)
    if TEMPDIR.exists():
        shutil.rmtree(TEMPDIR)
    UPLOAD_FOLDER = Path(TEMPDIR, "upload")
    UPLOAD_FOLDER.mkdir(parents=True, exist_ok=True)
    app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER

    files = request.files.getlist("sample")
    upload_files(files)
    PATH_SAMPLE = get_path_of_uploaded_files(UPLOAD_FOLDER, files)

    files = request.files.getlist("control")
    upload_files(files)
    PATH_CONTROL = get_path_of_uploaded_files(UPLOAD_FOLDER, files)

    files = request.files.getlist("allele")
    upload_files(files)
    PATH_ALLELE = get_path_of_uploaded_files(UPLOAD_FOLDER, files)

    genome = request.form.get("genome")
    threads = request.form.get("threads")
    if threads is None:
        threads = 1

    data = []
    for sample in PATH_SAMPLE:
        row = {"sample": sample, "name": name, "control": PATH_CONTROL[0], "allele": PATH_ALLELE[0]}
        if genome:
            row["genome"] = genome
        data.append(row)

    output_path = Path("TEMPDIR", "upload", "batch.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

    arguments = {}
    arguments["file"] = str(Path(TEMPDIR, "upload", "batch.csv"))
    arguments["threads"] = threads
    arguments["debug"] = False

    main.execute_batch_mode(arguments)

    return f"""
    name={name}
    sample={PATH_SAMPLE}
    control={PATH_CONTROL}
    allele={PATH_ALLELE}
    genome={genome}
    threads={threads}
    arguments={arguments["file"]}
    """


def open_browser(PORT):
    webbrowser.open_new(f"http://localhost:{PORT}/")


def execute():
    PORT = find_free_port()
    print(f"Assess 'http://localhost:{PORT}/' if a browser does not automatically open.")
    Timer(1, open_browser, [PORT]).start()
    with redirect_stderr(open(os.devnull, "w")):
        serve(app, host="0.0.0.0", port=PORT)
