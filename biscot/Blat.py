from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, OrderedDict

import logging
import os
import subprocess


def run_blat(
    contig_1,
    contig_1_start,
    contig_1_end,
    contig_1_orientation,
    contig_2,
    contig_2_start,
    contig_2_end,
    contig_2_orientation,
    contigs_sequence_dict,
):
    """
    Launches Blat, then sorts the psl file and extracts the best match
    
    :param contig_1: Name of the first contig
    :type contig_1: str
    :param contig_1_start: Starting position of the sequence to extract
    :type contig_1_start: int
    :param contig_1_end: End position of the sequence to extract
    :type contig_1_end: int
    :param contig_1_orientation: Strand orientation of the first contig
    :type contig_1_orientation: str
    :param contig_2: Name of the second contig
    :type contig_2: str
    :param contig_2_start: Starting position of the sequence to extract
    :type contig_2_start: int
    :param contig_2_end: End position of the sequence to extract
    :type contig_2_end: int
    :param contig_2_orientation: Strand orientation of the second contig
    :type contig_2_orientation: str
    :param contigs_sequence_dict: Dict containing the contigs fasta sequence
    :type contigs_sequence_dict: dict(str, str)
    """

    with open("Phase_2/tmp_1.fasta", "w") as out:
        if contig_1_orientation == "+":
            out.write(
                SeqRecord(
                    Seq(
                        contigs_sequence_dict[contig_1][contig_1_start:contig_1_end],
                        generic_dna,
                    ),
                    id=contig_1,
                    description="",
                ).format("fasta")
            )
        else:
            out.write(
                SeqRecord(
                    Seq(
                        contigs_sequence_dict[contig_1][contig_1_start:contig_1_end],
                        generic_dna,
                    ),
                    id=contig_1,
                    description="",
                )
                .reverse_complement()
                .format("fasta")
            )

    with open("Phase_2/tmp_2.fasta", "w") as out:
        if contig_2_orientation == "+":
            out.write(
                SeqRecord(
                    Seq(
                        contigs_sequence_dict[contig_2][contig_2_start:contig_2_end],
                        generic_dna,
                    ),
                    id=contig_2,
                    description="",
                ).format("fasta")
            )
        else:
            out.write(
                SeqRecord(
                    Seq(
                        contigs_sequence_dict[contig_2][contig_2_start:contig_2_end],
                        generic_dna,
                    ),
                    id=contig_2,
                    description="",
                )
                .reverse_complement()
                .format("fasta")
            )

    _ = subprocess.run(
        [
            "blat",
            "-noHead",
            "-minScore=3000",
            "Phase_2/tmp_1.fasta",
            "Phase_2/tmp_2.fasta",
            "Phase_2/tmp.psl",
        ],
        stdout=open("Phase_2/blat.out", "w"),
        stderr=open("Phase_2/blat.err", "w"),
    )

    _ = subprocess.run(
        ["pslSort", "dirs", "Phase_2/tmp.psl.sorted", "/tmp", "Phase_2/tmp.psl"],
        stdout=open("Phase_2/pslSort.out", "w"),
        stderr=open("Phase_2/pslSort.err", "w"),
    )

    _ = subprocess.run(
        [
            "pslReps",
            "-nohead",
            "-singleHit",
            "Phase_2/tmp.psl.sorted",
            "Phase_2/tmp.psl.besthit",
            "Phase_2/tmp.psr",
        ],
        stdout=open("Phase_2/pslReps.out", "w"),
        stderr=open("Phase_2/pslReps.err", "w"),
    )


def parse_blat():
    """
    Parses the pslReps output file to get various information about the best hit
    
    :return: (true if hit found - False otherwise, max score, reference size, contig_1 end, contig_2 end)
    :rtype: tuple(bool, int, int, int)
    """

    with open("Phase_2/tmp.psl.besthit") as blat:
        max_score = 0
        max_score_line = ""

        for line in blat:
            line = line.split("\t")

            if int(line[0]) > max_score:
                max_score = int(line[0])
                max_score_line = line

        if max_score > 0:
            return (
                True,
                max_score,
                int(max_score_line[14]),
                int(max_score_line[16]),
                int(max_score_line[12]),
            )

        return (False, None, None, None)


def get_agp_changes(contigs_sequence_dict):
    """
    Gets the best blat hit and writes what changes have to be made to the AGP file to incorporate the changes
    
    :param contigs_sequence_dict: Dict containing contigs FASTA sequences
    :type contigs_sequence_dict: dict(str: str)
    """

    agp_lines = open("Phase_1/phase_1.agp").readlines()

    for i, line in enumerate(agp_lines):
        (
            scaffold_name,
            scaffold_start,
            _,
            _,
            sequence_type,
            *other_fields,
        ) = line.rstrip().split("\t")
        gap_size = int(other_fields[0]) if sequence_type == "N" else 100000000000000000

        if gap_size < 1000:
            (
                previous_contig,
                previous_contig_start,
                previous_contig_end,
                previous_contig_orientation,
            ) = (agp_lines[i - 1].rstrip().split("\t")[5:9])

            next_contig, next_contig_start, next_contig_end, next_contig_orientation = (
                agp_lines[i + 1].rstrip().split("\t")[5:9]
            )

            if (
                len(contigs_sequence_dict[previous_contig]) > 50000
                and len(contigs_sequence_dict[next_contig]) > 50000
            ):
                run_blat(
                    previous_contig,
                    int(previous_contig_end) - 30000,
                    int(previous_contig_end),
                    previous_contig_orientation,
                    next_contig,
                    int(next_contig_start),
                    int(next_contig_start) + 30000,
                    next_contig_orientation,
                    contigs_sequence_dict,
                )

                if os.stat("Phase_2/tmp.psl.besthit").st_size == 0:
                    logging.debug(
                        f"NO MATCH between {previous_contig} and {next_contig}"
                    )
                    continue

                (
                    match_found,
                    max_score,
                    ref_size,
                    aln_ref_end,
                    aln_query_end,
                ) = parse_blat()

                if match_found:
                    logging.debug(
                        f"MATCH between {previous_contig} and {next_contig}, score: {max_score}"
                    )
                    with open("Phase_2/agp_changes.txt", "a") as changes_file:
                        changes_file.write(
                            "%s\t%s\t%s\t%s\t%s\n"
                            % (
                                scaffold_name,
                                scaffold_start,
                                aln_ref_end,
                                aln_query_end,
                                ref_size,
                            )
                        )


def write_new_agp(new_agp_lines):
    """
    Writes a new AGP file incorporating the Blat changes
    
    :param new_agp_lines: List containing the new AGP lines
    :type new_agp_lines: list(list(str))
    """

    with open("scaffolds.agp", "w") as out:
        last_scaffold = ""
        last_position = 0
        last_id = 0

        for line in new_agp_lines:
            if line[0] != last_scaffold:
                last_scaffold = ""
                last_position = 0
                last_id = 0

            if int(line[1]) != last_position + 1:
                # New end = end - (start - last_pos)
                line[2] = str(int(line[2]) - (int(line[1]) - last_position))

                # New start = last_pos + 1
                line[1] = str(last_position + 1)

            line[3] = str(int(last_id) + 1)

            last_position = int(line[2])
            last_scaffold = line[0]
            last_id = line[3]

            out.write("\t".join(line))


def mute_agp_file():
    """
    Loads an AGP file and a changes file and modifies the AGP file lines to incorporate the changes
    """

    agp_file = open("Phase_1/phase_1.agp")
    changes_file = open("Phase_2/agp_changes.txt")

    changes_dict = defaultdict(list)
    for line in changes_file:
        line = line.strip().split("\t")
        changes_dict[line[0]].append((line[1], line[2], line[3], line[4]))
    changes_file.close()

    agp_dict = OrderedDict()
    for line in agp_file:
        line = line.split("\t")
        try:
            agp_dict[line[0]].append(line)
        except:
            agp_dict[line[0]] = []
            agp_dict[line[0]].append(line)

    new_agp_lines = []
    # Mutating lines with a change
    for scaffold in agp_dict:
        for i, line in enumerate(agp_dict[scaffold]):
            if scaffold in changes_dict:
                change_found = False
                for change in changes_dict[scaffold]:
                    # Verify that the position is the same and that we are not examining an 'N' line
                    if (change[0] == line[1]) and (
                        agp_dict[scaffold][i + 1][6] != "scaffold"
                    ):
                        # Change start of next line on scaffold
                        agp_dict[scaffold][i + 1][1] = str(
                            int(agp_dict[scaffold][i + 1][1]) + int(change[2])
                        )

                        # Change start of next line on contig
                        agp_dict[scaffold][i + 1][6] = str(
                            int(agp_dict[scaffold][i + 1][6]) + int(change[2])
                        )

                        # Change end of previous line on scaffold
                        new_agp_lines[-1][2] = str(
                            int(agp_dict[scaffold][i - 1][2])
                            - (int(change[3]) - int(change[1]))
                        )

                        # Change end of previous line on contig
                        new_agp_lines[-1][7] = str(
                            int(agp_dict[scaffold][i - 1][7])
                            - (int(change[3]) - int(change[1]))
                        )

                        change_found = True
                        logging.debug(
                            f"Removed gap (size: {line[5]}) at position {line[1]} on scaffold {scaffold}"
                        )
                        break

                if not change_found:
                    new_agp_lines.append(agp_dict[scaffold][i])

            else:
                new_agp_lines.append(agp_dict[scaffold][i])
    write_new_agp(new_agp_lines)


def blat_phase(contigs_sequence_dict):
    """
    Executes Blat and modifies an AGP file based on the mappings found
    
    :param contigs_sequence_dict: Dict containing contigs FASTA sequences
    :type contigs_sequence_dict: dict(str: str)
    """

    try:
        os.mkdir("Phase_2")
    except:
        logging.debug("Phase_2 directory already created")

    with open("Phase_2/agp_changes.txt", "w"):
        pass

    get_agp_changes(contigs_sequence_dict)
    mute_agp_file()
