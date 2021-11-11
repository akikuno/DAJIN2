import shutil
import os
import glob
import filecmp

"""
- Check whether the Control is cached.
- Load the control if the control is cached.
- Save the newly-generated control if the control is not chached.
"""


def exists(control: str, cachedir: str) -> bool:
    with open(control, "rb") as f:
        tmp_file = [f.readline() for _ in range(1000)]
    with open(os.path.join(cachedir, "cache_header2.txt"), "w") as f:
        f.write(str(tmp_file))
    file1 = os.path.join(cachedir, "cache_header1.txt")
    file2 = os.path.join(cachedir, "cache_header2.txt")
    return filecmp.cmp(file1, file2)

    # control_name = os.path.basename(control)
    # f = os.path.join(cachedir, control_name)
    # if os.path.exists(f):
    #     if os.path.getsize(control) == os.path.getsize(f):
    #         return True
    # return False


def save_header(control: str, chachedir: str) -> None:
    with open(control, "rb") as f:
        tmp_file = [f.readline() for _ in range(1000)]
    with open(os.path.join(chachedir, "cache_header1.txt"), "w") as f:
        f.write(str(tmp_file))


def load(cachedir: str, tmpdir: str) -> None:
    suffix = os.path.basename(tmpdir)
    control_files = glob.glob(os.path.join(cachedir, "*" + suffix))
    for f in control_files:
        shutil.copy(f, tmpdir)


def save(control: str, tmpdir: str, cachedir: str) -> None:
    control_name = os.path.basename(control).split(".", 1)[0]
    control_files = glob.glob(os.path.join(tmpdir, control_name) + "*")
    for f in control_files:
        shutil.copy(f, cachedir)

