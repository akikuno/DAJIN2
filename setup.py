import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt") as requirements_file:
    install_requirements = requirements_file.read().splitlines()

setuptools.setup(
    name="DAJIN2",
    version="0.1.3",
    author="Akihiro Kuno",
    author_email="akuno@md.tsukuba.ac.jp",
    description="One-step genotyping tools for Nanopore amplicon sequencing",
    long_description=long_description,
    url="https://github.com/akikuno/DAJIN2",
    install_requires=install_requirements,
    packages=setuptools.find_packages(where="src",),
    package_dir={"": "src"},
    entry_points={"console_scripts": ["DAJIN2=DAJIN2.DAJIN2:main"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
    ],
    python_requires=">=3.7",
)
