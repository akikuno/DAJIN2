import setuptools
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

with open("requirements.txt") as requirements_file:
    install_requirements = requirements_file.read().splitlines()

setuptools.setup(
    name="DAJIN2",
    version="0.3.4",
    author="Akihiro Kuno",
    author_email="akuno@md.tsukuba.ac.jp",
    description="One-step genotyping tools for targeted long-read sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/akikuno/DAJIN2",
    install_requires=install_requirements,
    packages=setuptools.find_packages(
        where="src",
    ),
    package_dir={"": "src"},
    entry_points={"console_scripts": ["DAJIN2=DAJIN2.main:execute"]},
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
