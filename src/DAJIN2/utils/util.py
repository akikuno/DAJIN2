import sys
import os


def check_dependencies(dependencies: list) -> None:
    for dep in dependencies:
        if not shutil.which(dep):
            print(f"{dep} is not found. Please install it.")
            sys.exit(1)


def make_directories(maindir: str, subdirs: list) -> None:
    os.makedirs(maindir, exist_ok=True)
    for sub in subdirs:
        dir = os.path.join(maindir, sub)
        os.makedirs(dir, exist_ok=True)
