from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

import coloredlogs
import os
import logging
import sys


_orig_print = print


def print(*args, **kwargs):
    """Unbuffured print function"""

    _orig_print(*args, flush=True, **kwargs)


def check_path(file_paths):
    """
    Checks if the files given as input exist
    
    :param file_paths: List of file paths
    :type file_paths: list(str)
    :raises FileNotFoundError: If a file is not found
    """

    for file_path in file_paths:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_path} file not found.")


def setup_logging(debug):
    """
    Setups the logging streams
    
    :param debug: True to activate debug logs
    :type debug: bool
    """

    coloredlogs.DEFAULT_FIELD_STYLES["asctime"] = {"color": "cyan", "bright": True}
    coloredlogs.DEFAULT_FIELD_STYLES["filename"] = {"color": "yellow", "bright": True}
    coloredlogs.DEFAULT_FIELD_STYLES["lineno"] = {"color": "yellow", "bright": True}
    coloredlogs.DEFAULT_FIELD_STYLES["levelname"] = {"bold": True}

    coloredlogs.DEFAULT_LEVEL_STYLES["INFO"] = {"color": "green", "bright": True}
    coloredlogs.DEFAULT_LEVEL_STYLES["DEBUG"] = {"color": "magenta", "bright": True}

    level = ""

    if debug:
        logging.basicConfig(
            format="[%(asctime)s]\t[%(filename)9.9s - %(lineno)03d]\t[%(levelname)5.5s]\t%(message)s",
            level=logging.DEBUG,
            filename="biscot.log",
        )
        level = "DEBUG"
    else:
        logging.basicConfig(
            format="[%(asctime)s]\t[%(filename)9.9s - %(lineno)03d]\t[%(levelname)5.5s]\t%(message)s",
            level=logging.INFO,
            filename="biscot.log",
        )
        level = "INFO"

    log_formatter = logging.Formatter(
        "[%(asctime)s]\t[%(filename)9.9s - %(lineno)03d]\t[%(levelname)5.5s]\t%(message)s"
    )
    root_logger = logging.getLogger()
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    coloredlogs.install(
        level=level,
        logger=root_logger,
        fmt="[%(asctime)s]\t[%(filename)9.9s - %(lineno)03d]\t[%(levelname)5.5s]\t%(message)s",
    )


def load_contigs(contigs_sequence_dict, contigs_path):
    """
    Extracts contig sequences from a FASTA file
    
    :param contigs_sequence_dict: Dict that will contain contigs FASTA sequence
    :type contigs_sequence_dict: dict(str: str)
    :param contigs_path: Path to a contigs FASTA file
    :type contigs_path: str
    """

    logging.info("Loading contigs fasta file")
    for record in SeqIO.parse(open(contigs_path), "fasta"):
        contigs_sequence_dict[record.id] = str(record.seq)


def agp_to_fasta(contigs_sequence_dict, agp_path, output_file):
    """
    Parses an AGP file and thanks to a dict containing contigs sequence, transforms it into a scaffolds FASTA File
    
    :param contigs_sequence_dict: Dict containing contigs FASTA sequence
    :type contigs_sequence_dict: dict(str: str)
    :param agp_path: Path to an AGP file
    :type agp_path: str
    :param output_file: Path to an output FASTA file
    :type output_file: str
    """

    logging.info("Converting agp file to fasta")
    scaffolds_sequence_dict = defaultdict(lambda: "")
    with open(agp_path) as agp:
        for line in agp:
            line = line.rstrip("\n").split("\t")

            seq_type = line[4]

            if seq_type == "W":
                scaffold_name, contig_name, contig_start, contig_end, orientation = (
                    line[0],
                    line[5],
                    int(line[6]),
                    int(line[7]),
                    line[8],
                )

                if orientation == "+":
                    scaffolds_sequence_dict[scaffold_name] += contigs_sequence_dict[
                        contig_name
                    ][contig_start:contig_end]
                else:
                    scaffolds_sequence_dict[scaffold_name] += str(
                        Seq(
                            contigs_sequence_dict[contig_name][contig_start:contig_end],
                            generic_dna,
                        ).reverse_complement()
                    )

            elif seq_type == "N":
                scaffold_name, gap_size = line[0], int(line[5])
                scaffolds_sequence_dict[scaffold_name] += "N" * gap_size

    sorted_scaffolds = sorted(
        scaffolds_sequence_dict.items(), key=lambda d: len(d[1]), reverse=True
    )
    with open(output_file, "w") as out:
        for scaffold, sequence in sorted_scaffolds:
            record = SeqRecord(Seq(sequence, generic_dna), id=scaffold, description='')
            out.write(record.format("fasta"))
