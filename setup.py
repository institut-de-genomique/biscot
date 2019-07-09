import setuptools
from pathlib import Path
from typing import List
from distutils.command.clean import clean


def parse_requirements(filename: str) -> List[str]:
    """Return requirements from requirements file."""
    # Ref: https://stackoverflow.com/a/42033122/
    requirements = (Path(__file__).parent / filename).read_text().strip().split('\n')
    requirements = [r.strip() for r in requirements]
    requirements = [r for r in sorted(requirements) if r and not r.startswith('#')]
    return requirements


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='biscot',  
    version='0.1',
    scripts=['biscot.py'] ,
    author="Benjamin Istace, Caroline Belser, Jean-Marc Aury",
    author_email="bistace@genoscope.cns.fr",
    description="Bionano SCaffolding Correction Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/institut-de-genomique/biscot",
    packages=["biscot_libs"],
    install_requires=parse_requirements('requirements.txt'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CeCILL",
        "Operating System :: Linux",
    ]
)
