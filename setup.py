from setuptools import setup, find_packages
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="biscot",
    packages=["biscot"],
    version="2.3",
    entry_points={
        "console_scripts": [
            "biscot = biscot.biscot:run",
        ],
    },
    license="CeCILL",
    description="Bionano SCaffolding Correction Tool",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Benjamin Istace, Caroline Belser, Jean-Marc Aury",
    author_email="bistace@genoscope.cns.fr",
    url="https://github.com/institut-de-genomique/biscot",
    download_url="https://github.com/institut-de-genomique/biscot/archive/v2.3.tar.gz",
    keywords=["Genome", "Scaffolding", "Bionano",],
    install_requires=["argparse", "biopython", "coloredlogs",],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Programming Language :: Python :: 3.6",
    ],
)
