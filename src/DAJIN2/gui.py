import socket
import webbrowser
from contextlib import closing
from pathlib import Path
from threading import Timer

import pandas as pd
from flask import Flask, render_template, request
from waitress import serve
from werkzeug.utils import secure_filename

from . import batch


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
    UPLOAD_FOLDER = Path("DAJINResults", ".tempdir", name, "upload")
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

    df = pd.DataFrame({"sample": PATH_SAMPLE})
    df["name"] = name
    df["control"] = PATH_CONTROL[0]
    df["allele"] = PATH_ALLELE[0]
    if genome:
        df["genome"] = genome
    df.to_csv(Path("DAJINResults", ".tempdir", name, "upload", "batch.csv"), index=False)

    arguments = dict()
    arguments["file"] = str(Path("DAJINResults", ".tempdir", name, "upload", "batch.csv"))
    arguments["threads"] = threads
    arguments["debug"] = False
    batch.execute(arguments)

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
    webbrowser.open_new(f"http://127.0.0.1:{PORT}/")


def execute():
    PORT = find_free_port()
    print(f"Assess 'http://127.0.0.1:{PORT}/' if a browser does not automatically open.")
    Timer(1, open_browser, [PORT]).start()
    serve(app, host="0.0.0.0", port=PORT)
