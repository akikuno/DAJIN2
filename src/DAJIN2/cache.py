import shutil
import os
import glob


def exists(control: str, cachedir: str) -> bool:
    control_name = os.path.basename(control)
    f = os.path.join(cachedir, control_name)
    return os.path.exists(f)


def is_same_filesize(control: str, cachedir: str) -> bool:
    control_name = os.path.basename(control)
    f = os.path.join(cachedir, control_name)
    return os.path.getsize(control) == os.path.getsize(f)


def save(control: str, tmpdir: str, cachedir: str) -> None:
    control_name = os.path.basename(control).split(".", 1)[0]
    control_files = glob.glob(os.path.join(tmpdir, control_name) + "*")
    for f in control_files:
        shutil.copy(f, cachedir)


def load(cachedir: str, tmpdir: str) -> None:
    suffix = os.path.basename(tmpdir)
    control_files = glob.glob(os.path.join(cachedir, "*" + suffix))
    for f in control_files:
        shutil.copy(f, tmpdir)
