from __future__ import annotations
from urllib.request import urlopen
from urllib.error import URLError


def check_url(urls: list[str]) -> tuple[str, bool]:
    flag_fail = False
    url = ""
    for url in urls:
        try:
            _ = urlopen(url)
        except URLError:
            flag_fail = True
        else:
            flag_fail = False
            break
    return url, flag_fail


def check_genome(genome: str, ucsc_url: str):
    url = f"{ucsc_url}/cgi-bin/das/{genome}/dna?segment=1:1,10"
    _, flag_fail = check_url([url])
    if flag_fail:
        raise AttributeError(f"{genome} is not listed in UCSC genome browser")
