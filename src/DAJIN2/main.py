import argparser
import shutil
import sys


def check_deps(dependencies):
    for dep in dependencies:
        if not shutil.which(dep):
            print(f"{dep} is not found. Please install it.")
            sys.exit(1)


def main():
    check_deps(["minimap2", "samtools"])
    sample, control, output, genome, debug, threads = argparser.parse()

    if debug is False:
        shutil.rmtree(".tmpDAJIN")


if __name__ == "__main__":
    main()
