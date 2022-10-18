import os
import re
from pathlib import Path
import http.server
import socketserver
import webbrowser
import socket
from contextlib import closing
from jinja2 import Template, Environment, FileSystemLoader


def find_free_port():
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


def execute(name: str):
    DIR_IGVJS = Path("DAJINResults", name, ".igvjs")
    if not DIR_IGVJS.exists():
        raise FileNotFoundError("BAM files for DAJIN view is not found. Execute DAJIN first.")
    env = Environment(loader=FileSystemLoader("./", encoding="utf8"))
    template = env.get_template("src/DAJIN2/template_igvjs.html")
    params_genome = {"genome": {"exist": False}}
    params_reference = dict()
    path_genome = Path(DIR_IGVJS, "genome_symbol.txt")
    if path_genome.exists():
        path_coodinates = Path(DIR_IGVJS, "genome_coodinates.jsonl")
        GENOME = path_genome.read_text().strip()
        CHROME, START, END, _ = eval(path_coodinates.read_text().strip()).values()
        params_genome = {"genome": {"exist": True, "genome": GENOME, "locus": f"{CHROME}:{START}-{END}"}}
    else:
        params_reference = {"reference": {"urlfa": "control.fasta", "urlfai": "control.fasta" + ".fai"}}

    bamnames = []
    bamurls = []
    baiurls = []
    for bam in Path("DAJINResults", name, ".igvjs").iterdir():
        if not re.search(r"bam$", str(bam)):
            continue
        bamnames.append(bam.stem)
        bamurls.append(str(bam.name))
        baiurls.append(str(bam.name) + ".bai")
    bamnames.sort()
    bamurls.sort()
    baiurls.sort()
    contents = [{"samplename": n, "urlbam": b, "urlbai": i} for n, b, i in zip(bamnames, bamurls, baiurls)]

    params_tracks = """{
        "tracks":
            {{ contents }},
    }"""
    params_tracks = Template(params_tracks)
    params_tracks = eval(params_tracks.render(contents=contents))

    params_genome.update(**params_reference, **params_tracks)
    params = params_genome.copy()

    HTML_IGVJS = template.render(params)
    Path(DIR_IGVJS, "index.html").write_text(HTML_IGVJS)
    PORT = find_free_port()
    Handler = http.server.SimpleHTTPRequestHandler
    os.chdir(Path("DAJINResults", name, ".igvjs"))
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        print(f"serving at port: http://127.0.0.1:{PORT}")
        httpd.serve_forever()
        webbrowser.open(f"http://127.0.0.1:{PORT}", autoraise=True)
    os.chdir("../../../")
