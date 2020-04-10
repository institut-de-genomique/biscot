# Bionano Scaffolding Correction Tool (BiSCoT)
BiSCoT is a tool that aims to improve the contiguity of scaffolds and contigs generated after a Bionano scaffolding. It looks for enzymatic labelling sites on contigs. If two distinct contigs share labels, BiSCoT merges them at the last shared site.

Biorxiv preprint : [link](https://www.biorxiv.org/content/10.1101/674721v1 "BiSCoT Biorxiv preprint")

For additional documentation and test cases, please have a look at the [wiki](https://github.com/institut-de-genomique/biscot/wiki "Wiki link").

In case of troubles when using or installing the software, please open up an issue by clicking [here](https://github.com/institut-de-genomique/biscot/issues/new "Github issue page").

## Installation

BiSCoT comes in the form of a Python3 script with some Python and software dependencies. In order to run it correctly, you will need :
- Python 3 (tested with Python 3.6)
- the Biopython and Argparse python modules (both installed automatically with BiSCoT)
- the BLAT aligner if you plan on using the aggressive scaffolding mode that is based on shared labels and sequence similarity (BiSCoT was tested with the [v36](https://hgwdev.gi.ucsc.edu/~kent/src/blatSrc36.zip "BLAT v36") version of BLAT)

BiSCoT is available on [PyPI](https://pypi.org/ "PyPI") and can be installed with the following command:
```
pip install biscot
```


## Usage

BiSCoT was designed to improve a prior Bionano scaffolding so it needs a few files generated during this step :
- one CMAP file (--cmap-ref argument) describing the positions of enzymatic labelling sites on the reference genome (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_NGS_contigs_HYBRID_Export_r.cmap` in the case of a double hybrid scaffolding)
- one CMAP file per enzyme (--cmap-1 and --cmap-2 arguments) describing the positions of enzymatic labelling sites on the contigs (filenames usually look like this : `E_CTTAAG_Q_NGScontigs_A_HYBRID_q.cmap` for DLE1 and `E_GCTCTTC_Q_NGScontigs_A_HYBRID_q.cmap` for BspQI)
- a KEY file (--key argument) describing the names of the contig maps related to their FASTA file header names (filename usually looks like this `*.cut_CTTAAG_GCTCTTC_0kb_0labels_key.txt`)
- one XMAP file per enzyme (--xmap-1 and --xmap-2 arguments) describing the alignments of contig labels on the anchor (filename usually looks like this `E_CTTAAG_Q_NGScontigs_A_HYBRID.xmap` for DLE1 and `E_GCTCTTC_Q_NGScontigs_A_HYBRID.xmap` for BspQI)
- the contigs FASTA file (--contigs argument) that was used for the scaffolding

A typical execution of BiSCoT should look like this :
```
# If the previous command returned an error, modify your $PATH to include the blat directory
export PATH=Path/To/Blat:$PATH

# Execute BiSCoT
biscot.py --cmap-ref cmap_reference.cmap \\
    --cmap-1 cmap_dle.cmap \\
    --cmap-2 cmap_bspq1.cmap \\
    --xmap-1 xmap_dle.xmap \\
    --xmap-2 xmap_bspqi.xmap \\
    --key key.txt \\
    --contigs contigs.fasta \\
    --output biscot
```

If everything went fine, a `biscot` directory should have been created. Inside, you will find two output files :
- a `scaffolds.fasta` containing the new scaffolds 
- a `scaffolds.agp` file containing the changes made to contigs

If you would like to change the name/path of the output directory, you can do so with the `--output` argument.

### Optional arguments
- `--xmap-2enz` argument is used to provide the final XMAP file containing the mappings of labels of both enzymes. This argument is useful (and recommended) to ensure that no mapping has been missed inside one of the individual XMAP file. Usually, this file's name looks like this: `*.cut_CTTAAG_GCTCTTC_0kb_0labels_NGS_contigs_HYBRID_Export.xmap_sorted.xmap`.
- `--only-confirmed-pos` argument is used so that only mappings contained in the `--xmap-2enz` file are retained. Indeed, by using one XMAP file per enzyme, contigs can be placed two times in the final assembly. This ensures that contigs are only placed one time and that created scaffolds are validated by both enzymes.
- `--aggressive` enables the sequence similarity scaffolding. In a first phase, BiSCoT will search similarities between contigs based on label mappings. If this parameter is set, BiSCoT will search for sequence similarity to close gaps created by the first step.
