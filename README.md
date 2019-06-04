# Bionano Scaffolding Correction Tool (BiSCoT)
BiSCoT is a tool that aims to improve the contiguity of scaffolds and contigs generated after a Bionano scaffolding. It looks for enzymatic labelling sites on contigs. If two distinct contigs share labels, BiSCoT merges them at the last shared site.

Biorxiv preprint : [XXX](https://www.biorxiv.org/ "BiSCoT Biorxiv preprint")

In case of troubles when using or installing the software, please open up an issue by clicking [here](https://github.com/institut-de-genomique/biscot/issues/new "Github issue page").

## Installation

BiSCoT comes in the form of a Python3 script with some Python and software dependencies. In order to run it correctly, you will need :
- Python 3 (tested with Python 3.6)
- virtualenv (bundled with Python, starting from Python 3.5)
- the Biopython and Argparse python modules (both installed automatically with BiSCoT)
- the BLAT aligner (BiSCoT was tested with the [v36](https://hgwdev.gi.ucsc.edu/~kent/src/blatSrc36.zip "BLAT v36") version of BLAT)

First, clone this github repository :<br>
`git clone https://github.com/institut-de-genomique/biscot.git`<br>

Then, move to the downloaded repository, create a virtualenv and install BiSCoT and its Python dependencies :<br>
```
cd biscot
virtualenv venv
source venv/bin/activate
python setup.py install
python setup.py clean
```
When you're done with using BiSCoT and you want to leave the virtualenv, simply type :
```
deactivate
```

## Usage

BiSCoT was designed to improve a prior Bionano scaffolding so it needs a few files generated during this step :
- a CMAP file of the reference genome (--cmap-ref argument) describing the positions of enzymatic labelling sites on the anchor (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_NGS_contigs_HYBRID_Export_r.cmap` in the case of a double hybrid scaffolding)
- one or two CMAP files of the contigs fasta file (--cmap-dle and --cmap-bspq1 arguments) describing the positions of enzymatic labelling sites on the contigs (filenames usually look like this : `E_CTTAAG_Q_NGScontigs_A_HYBRID_q.cmap` for DLE1 and `E_GCTCTTC_Q_NGScontigs_A_HYBRID_q.cmap` for BspQI)
- a KEY file (--key argument) describing the names of the contig maps related to their FASTA file header names (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_key.txt`)
- a XMAP file (--xmap argument) describing the alignments of contig labels on the anchor (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_NGS_contigs_HYBRID_Export.xmap_sorted.xmap`)
- the contigs FASTA file (--contigs argument) that was used for the scaffolding

Final command line should look like this :
```
# Activate the virtualenv, replace BISCOT_DIR with the path to your BiSCoT installation
source BISCOT_DIR/venv/bin/activate

# Make sure that BLAT is in your $PATH
which blat

# If the previous command returned an error, modify your $PATH to include the blat directory
export PATH=Path/To/Blat:$PATH

# Execute BiSCoT
biscot.py --cmap-ref cmap_reference.cmap \\
    --cmap-1 cmap_dle.cmap \\
    --cmap-2 cmap_bspq1.cmap \\
    --xmap xmap.xmap \\
    --key key.txt \\
    --contigs contigs.fasta \\
    --output biscot
    
# Exit the virtualenv
deactivate
```

If everything went fine, a `biscot` directory should have been created. Inside, you will find two output files :
- a `scaffolds.fasta` containing the new scaffolds 
- a `scaffolds.agp` file containing the changes made to contigs. 
If you would like to change the name/path of the output directory, you can do so with the `--output` argument.
