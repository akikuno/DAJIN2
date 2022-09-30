from setuptools import setup, find_packages

with open("requirements.txt") as requirements_file:
    install_requirements = requirements_file.read().splitlines()

from setuptools import setup

setup(
    name="DAJIN2",
    version="1.0.0",
    description="One-step genotyping tools for Nanopore amplicon sequencing",
    author="Akihiro Kuno",
    author_email="akuno@md.tsukuba.ac.jp",
    url="https://github.com/akikuno/DAJIN2",
    packages=find_packages(),
    install_requires=install_requirements,
    entry_points={"console_scripts": ["DAJIN2=src.DAJIN2.DAJIN2:main"]},
    # package_dir={"DAJIN2": "src/DAJIN2",},
    # package_data={"DAJIN2": ["core/*.py",]},
)
