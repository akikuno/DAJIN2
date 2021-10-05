import shutil
import os
import glob


def save_fastq(control: str, cachedir: str) -> None:
    shutil.copy(control, cachedir)


def exists(control: str, cachedir: str) -> bool:
    control_name = os.path.basename(control)
    f = os.path.join(cachedir, control_name)
    return os.path.exists(f)


def is_same(control: str, cachedir: str) -> bool:
    control_name = os.path.basename(control)
    f = os.path.join(cachedir, control_name)
    os.path.getsize(control) == os.path.getsize(f)
    shutil.copy(p, dst_path)


def save(control: str, cachedir: str) -> None:
    control_name = os.path.basename(control).split(".", 1)[0]
    control_files = glob.glob(os.path.join(tmpdir, control_name) + "*")
    for f in control_files:
        shutil.copy(f, cachedir)


def load(control: str, cachedir: str) -> None:
    pass
