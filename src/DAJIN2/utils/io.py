import gzip


def fread(filepath: str) -> list:
    """
    read gzip or text file
    """
    with gzip.open(filepath, "rt") as f:
        try:
            out = f.read().splitlines()
        except:
            with open(filepath) as f:
                out = f.read().splitlines()
    return out


def fwrite(infile: list, outfilepath: str) -> None:
    """
    save gzipped file
    """
    with gzip.open(outfilepath, "wt") as f:
        for i in infile:
            f.write("\n".join(i) + "\n")

