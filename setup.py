from distutils.core import setup

setup(
    name="biscot",
    packages=["biscot"],
    version="2.0",
    license="CeCILL",
    description="Bionano SCaffolding Correction Tool",
    author="Benjamin Istace, Caroline Belser, Jean-Marc Aury",
    author_email="bistace@genoscope.cns.fr",
    url="https://github.com/institut-de-genomique/biscot",
    download_url="https://github.com/institut-de-genomique/biscot/archive/v2.0.tar.gz",
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
