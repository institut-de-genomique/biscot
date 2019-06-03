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

## Usage

BiSCoT was designed to improve a prior Bionano scaffolding so it needs a few files generated during this step :
- one CMAP file of the reference genome (argument --cmap-ref) describing the positions of enzymatic labelling sites on the anchor (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_NGS_contigs_HYBRID_Export_r.cmap` in the case of a double hybrid scaffolding)
- one or two CMAP files of the contigs fasta file (arguments --cmap-dle and --cmap-bspq1) describing the positions of enzymatic labelling sites on the contigs (filenames usually look like this : `E_CTTAAG_Q_NGScontigs_A_HYBRID_q.cmap` for DLE1 and `E_GCTCTTC_Q_NGScontigs_A_HYBRID_q.cmap` for BspQI. If you don't find these two files in the directory you exported from the Bionano Access, see the section 'Finding contigs CMAP files' of this README)
- one KEY file (argument --key) describing the names of the contig maps related to their FASTA file header names (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_key.txt`)
- one XMAP file (argument --xmap) describing the alignments of contig labels on the anchor (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_NGS_contigs_HYBRID_Export.xmap_sorted.xmap`)
- the contigs FASTA file (argument --contigs) that was used for the scaffolding

Final command line should look like this :
```
source BISCOT_DIR/venv/bin/activate
biscot.py --cmap-ref cmap_refeference.cmap \\
    --cmap-bspq1 cmap_bspq1.cmap \\
    --cmap-dle cmap_dle.cmap
    --xmap xmap.xmap
    --key key.txt
    --contigs contigs.fasta
deactivate
```

    
