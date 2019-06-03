# Bionano SCaffolding Correction Tool (BiSCoT)
BiSCoT is a tool that aims to improve the contiguity of scaffolds and contigs generated after a Bionano scaffolding. In order to do that, it looks for enzymatic labelling sites on contigs and merges them if they share labels.

Biorxiv preprint : [XXX](https://www.biorxiv.org/ "BiSCoT Biorxiv preprint")

In case of troubles when using/installing the software, please open up an issue by clicking [here](https://github.com/institut-de-genomique/biscot/issues/new "Github issue page").

## Installation

BiSCoT comes in the form of a Python3 script with some (not much) Python and software dependencies. In order to run it correctly, you will need :
- Python 3 (BiSCoT was tested with Python 3.6)
- the Biopython and Argparse python modules (both installed automatically with BiSCoT)
- the BLAT aligner (BiSCoT was tested with the [v36](https://hgwdev.gi.ucsc.edu/~kent/src/blatSrc36.zip "BLAT v36") version of BLAT)

First, clone this github repository :<br>
`git clone https://github.com/institut-de-genomique/biscot.git`<br>

Then, move to the downloaded repository, create a virtualenv (included with every python distribution, starting from Python 3.5) and install BiSCoT and its Python dependencies :<br>
```
cd biscot
virtualenv venv
source venv/bin/activate
python setup.py install
python setup.py clean
```
When you're done with using BiSCoT, simply type `deactivate` to exit the virtualenv.
