from biscot import Key
from biscot import Map

from collections import defaultdict
import copy
import logging
import math


class Alignment:
    def __init__(
        self,
        map_id,
        reference_id,
        map_start,
        map_end,
        reference_start,
        reference_end,
        channel,
        orientation,
        alignments,
    ):
        self.map_id = map_id
        self.map_start = map_start
        self.map_end = map_end

        self.reference_id = reference_id
        self.reference_start = reference_start
        self.reference_end = reference_end

        self.channel = channel
        self.orientation = orientation
        self.alignment = copy.deepcopy(alignments)

    def get_corresponding_contig_map_label(self, reference_label_id):
        """
        Extracts a contig map label that was aligned against specified anchor label
        
        :param reference_label_id: Anchor map label id
        :type reference_label_id: int
        :return: Contig map label id
        :rtype: int
        """

        for reference, contig, _ in self.alignment:
            if reference == reference_label_id:
                return contig

    def update_alignments(self, contig_maps_dict):
        """
        Parses contig maps and removes alignment when aln.map_start < self.map_start 
        or aln.map_end > self.map_end
        
        :param contig_maps_dict: Dict containing contig Map objects
        :type contig_maps_dict: dict(int, Map)
        """

        alignments_to_remove = set()
        for reference, contig, channel in self.alignment:
            position = contig_maps_dict[self.map_id].get_label_position(contig, channel)
            if (position < min(self.map_start, self.map_end)) or (
                position > max(self.map_end, self.map_start)
            ):
                alignments_to_remove.add((reference, contig, channel))

        self.alignment = [i for i in self.alignment if i not in alignments_to_remove]

    def __str__(self):
        """
        Text representation of an Alignment object
        
        :return: String representing the Alignment object
        :rtype: str
        """

        txt = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
            self.reference_id,
            self.reference_start,
            self.reference_end,
            self.map_id,
            self.map_start,
            self.map_end,
            self.orientation,
            self.channel,
        )
        return txt


def line_to_alignment(line, channel):
    """
    Converts an xmap line to an Alignment object
    
    :param line: A line of an xmap file
    :type line: str
    :param channel: Enzyme channel to consider
    :type channel: int
    :return: An alignment object
    :rtype: Alignment
    """
    line = line.rstrip("\n").split("\t")

    map_id = int(line[1])
    map_start = int(line[3].split(".")[0])
    map_end = int(line[4].split(".")[0])

    reference_id = int(line[2])
    reference_start = int(line[5].split(".")[0])
    reference_end = int(line[6].split(".")[0])

    orientation = line[7]
    alignment = line[13]

    alignment = alignment.replace(")", "")
    alignment_list = []
    for label_couples in alignment.split("(")[1:]:
        label_couples = label_couples.split(",")
        reference = int(label_couples[0])
        contig = int(label_couples[1])
        alignment_list.append((reference, contig, 1))

    aln = Alignment(
        map_id,
        reference_id,
        map_start,
        map_end,
        reference_start,
        reference_end,
        channel,
        orientation,
        alignment_list,
    )

    return aln


def parse_xmap(
    reference_maps_dict,
    xmap_1_path,
    xmap_2_path,
    deleted_xmap_records,
    xmap_two_enzymes_path="",
    only_confirmed_positions=False,
):
    """
    Parses from one to three xmaps and converts lines to Alignment objects
    
    :param reference_maps_dict: Dict containing anchor Map obecjts
    :type reference_maps_dict: dict(int, Map)
    :param xmap_1_path: Path to the first xmap_file
    :type xmap_1_path: str
    :param xmap_2_path: Path to the second xmap file
    :type xmap_2_path: str
    :param deleted_xmap_records: Dict containing Alignment objects that were deleted due to a larger alignment being found
    :type deleted_xmap_records: dict(int, Alignment)
    :param xmap_two_enzymes_path: Path to the 2-enzyme xmap file, defaults to ""
    :type xmap_two_enzymes_path: str, optional
    :param only_confirmed_positions: If True, only alignments contained in xmap_1 or xmap_2 AND in xmap_2enzymes will be conserved, defaults to False
    :type only_confirmed_positions: bool, optional
    """

    xmap_1_file = open(xmap_1_path)

    for line in xmap_1_file:
        if not line.startswith("#"):
            aln = line_to_alignment(line, 1)
            reference_maps_dict[aln.reference_id].add_alignment(aln)
    xmap_1_file.close()

    if xmap_2_path:
        xmap_2_file = open(xmap_2_path)

        for line in xmap_2_file:
            if not line.startswith("#"):
                aln = line_to_alignment(line, 2)

                if aln.map_id in reference_maps_dict[aln.reference_id].contig_maps:
                    aln_1 = None
                    for alignment in reference_maps_dict[aln.reference_id].alignments:
                        if aln.map_id == alignment.map_id:
                            aln_1 = alignment

                    size_aln_1 = math.fabs(aln_1.map_start - aln_1.map_end)
                    size_aln_2 = math.fabs(aln.map_start - aln.map_end)
                    if size_aln_2 > size_aln_1:
                        deleted_xmap_records[aln.map_id] = copy.deepcopy(aln_1)
                        for i, alignment in enumerate(
                            reference_maps_dict[aln.reference_id].alignments
                        ):
                            if aln.map_id == alignment.map_id:
                                reference_maps_dict[aln.reference_id].alignments[
                                    i
                                ] = copy.deepcopy(aln)
                                break
                    else:
                        deleted_xmap_records[aln.map_id] = copy.deepcopy(aln)
                else:
                    reference_maps_dict[aln.reference_id].add_alignment(
                        copy.deepcopy(aln)
                    )
        xmap_2_file.close()

    if xmap_two_enzymes_path:
        xmap_two_enzymes_file = open(xmap_two_enzymes_path)
        scaffold_names = set()
        verified_contigs = defaultdict(list)

        for line in xmap_two_enzymes_file:
            if not line.startswith("#"):
                aln = line_to_alignment(line, 1)

                scaffold_names.add(aln.reference_id)

                for i in range(len(reference_maps_dict[aln.reference_id].alignments)):
                    if (
                        reference_maps_dict[aln.reference_id].alignments[i].map_id
                        == aln.map_id
                    ):
                        reference_maps_dict[aln.reference_id].alignments[
                            i
                        ].reference_start = aln.reference_start
                        reference_maps_dict[aln.reference_id].alignments[
                            i
                        ].reference_end = aln.reference_end
                        reference_maps_dict[aln.reference_id].alignments[
                            i
                        ].map_start = aln.map_start
                        reference_maps_dict[aln.reference_id].alignments[
                            i
                        ].map_end = aln.map_end
                        break

                if aln.map_id not in reference_maps_dict[aln.reference_id].contig_maps:
                    reference_maps_dict[aln.reference_id].add_alignment(aln)

                verified_contigs[aln.reference_id].append(aln.map_id)

        unverified_scaffolds = []
        for scaffold_id in reference_maps_dict:
            if scaffold_id not in scaffold_names:
                unverified_scaffolds.append(scaffold_id)
        for scaffold_id in unverified_scaffolds:
            reference_maps_dict.pop(scaffold_id)

        if only_confirmed_positions:
            logging.info("Removing unverified alignments")
            alignments_to_pop = []
            for scaffold in reference_maps_dict:
                for alignment in reference_maps_dict[scaffold].alignments:
                    if alignment.map_id not in verified_contigs[scaffold]:
                        alignments_to_pop.append(alignment)
                        logging.debug(
                            f"Unverified position - Alignment of map {alignment.map_id} on anchor {alignment.reference_id} will be removed"
                        )

            for scaffold in reference_maps_dict:
                reference_maps_dict[scaffold].alignments = [
                    alignment
                    for alignment in reference_maps_dict[scaffold].alignments
                    if alignment not in alignments_to_pop
                ]

        xmap_two_enzymes_file.close()


def get_shared_labels(aln_1, aln_2):
    """
    Parses two Alignments objects and returns anchor map label ids for which both contig maps are aligned to
    
    :param aln_1: An alignment object
    :type aln_1: Alignment
    :param aln_2: An Alignment object
    :type aln_2: Alignment
    :return: A tuple of lists that contains reference label ids which corresponds to the overlap between aln_1 and aln_2. One list for each channel.
    :rtype: tuple(list(int), list(int))
    """

    labels_aln_1 = {1: set(), 2: set()}
    labels_aln_2 = {1: set(), 2: set()}

    for reference, _, channel in aln_1.alignment:
        labels_aln_1[channel].add(reference)

    for reference, _, channel in aln_2.alignment:
        labels_aln_2[channel].add(reference)

    intersection_channel_1 = labels_aln_1[1].intersection(labels_aln_2[1])
    intersection_channel_2 = labels_aln_1[2].intersection(labels_aln_2[2])

    return (sorted(intersection_channel_1), sorted(intersection_channel_2))


def get_leftmost_label(label_list, channel, reference_map):
    """
    Extracts all label positions from a label list ids and reference maps and returns the one that has the minimum position on the anchor
    
    :param label_list: List of label ids
    :type label_list: list(int)
    :param channel: Enzyme channel to consider to extract label position
    :type channel: int (1 or 2)
    :param reference_map: Dict containing anchor Map objects
    :type reference_map: dict(integer: Map)
    :return: Returns the label id that satisfies label_position = min(all_label_positions)
    :rtype: int
    """

    leftmost_label = None
    min_pos = 10000000000

    for label in label_list:
        try:
            pos = reference_map.get_label_position(label, channel)
        except:
            logging.debug(
                f"Didn't find label {label} channel {channel} on reference map {reference_map.map_id}"
            )
            break
        if pos < min_pos:
            leftmost_label = (label, pos)
            min_pos = pos

    return leftmost_label


def get_rightmost_label(label_list, channel, reference_map):
    """
    Extracts all label positions from a label list ids and reference maps and returns the one that has the maximum position on the anchor
    
    :param label_list: List of label ids
    :type label_list: list(int)
    :param channel: Enzyme channel to consider to extract label position
    :type channel: int (1 or 2)
    :param reference_map: Dict containing anchor Map objects
    :type reference_map: dict(int, Map)
    :return: Returns the label id that satisfies label_position = max(all_label_positions)
    :rtype: int
    """

    rightmost_label = None
    max_pos = 0

    for label in label_list:
        try:
            pos = reference_map.get_label_position(label, channel)
        except:
            logging.debug(
                f"Didn't find label {label} channel {channel} on reference map {reference_map.map_id}"
            )
            break
        if pos > max_pos:
            rightmost_label = (label, pos)
            max_pos = pos

    return rightmost_label


def solve_containment(aln_couple, reference_maps_dict, contig_maps_dict, key_dict):
    """
    |  Tries to integrate a small map into a larger one.
    |  Let's consider a Map 1 that is aligned on the reference from position 1 to 100 and a Map 2 that is aligned on the reference from position 25 to 75.
    |  The goal of this function is to break alignment of Map 1 into two alignments (1-25 and 75-100).
    
    :param aln_couple: Two Alignment objects. The first one being the 'small alignment' and the second, the 'large alignment'
    :type aln_couple: tuple(Alignment, Alignment)
    :param reference_maps_dict: Dict of anchor Map objects
    :type reference_maps_dict: dict(int, Map)
    :param contig_maps_dict: Dict of contig Map objects
    :type contig_maps_dict: dict(int, Map)
    :param key_dict: Dict containing the correspondance between Map objects and actual sequences
    :type key_dict: dict((int, int, int), (str, int, int, int))
    """

    large_aln = aln_couple[0]
    short_aln = aln_couple[1]

    if (short_aln.reference_start - 10000 < large_aln.reference_start) or (
        short_aln.reference_end + 10000 > large_aln.reference_end
    ):
        logging.debug(
            f"Map {short_aln.map_id} is contained at the extremities of map {large_aln.map_id} (Anchor {short_aln.reference_id}), removing alignment of map {short_aln.map_id}"
        )

        for i in range(0, len(reference_maps_dict[short_aln.reference_id].alignments)):
            if reference_maps_dict[short_aln.reference_id].alignments[i] == short_aln:
                reference_maps_dict[short_aln.reference_id].alignments.pop(i)
                break

    else:
        intersection = get_shared_labels(short_aln, large_aln)
        if len(intersection[0]) + len(intersection[1]) <= 2:
            logging.debug(
                f"Not enough shared labels to integrate map {short_aln.map_id} in map {large_aln.map_id}, removing alignment of map {short_aln.map_id}"
            )

            for i in range(
                0, len(reference_maps_dict[short_aln.reference_id].alignments)
            ):
                if (
                    reference_maps_dict[short_aln.reference_id].alignments[i]
                    == short_aln
                ):
                    reference_maps_dict[short_aln.reference_id].alignments.pop(i)
                    break

        else:
            label_1_start = get_leftmost_label(
                intersection[0], 1, reference_maps_dict[large_aln.reference_id]
            )
            label_2_start = get_leftmost_label(
                intersection[1], 2, reference_maps_dict[large_aln.reference_id]
            )
            label_1_end = get_rightmost_label(
                intersection[0], 1, reference_maps_dict[large_aln.reference_id]
            )
            label_2_end = get_rightmost_label(
                intersection[1], 2, reference_maps_dict[large_aln.reference_id]
            )

            if (
                not label_1_start
                and not label_1_end
                and not label_2_start
                and not label_2_end
            ):
                logging.debug(
                    f"Couldn't find label boundaries for alignment of map {short_aln.map_id} contained in map {large_aln.map_id} on anchor {short_aln.reference_id}, removing alignment of map {short_aln.map_id}"
                )

                for i in range(
                    0, len(reference_maps_dict[short_aln.reference_id].alignments)
                ):
                    if (
                        reference_maps_dict[short_aln.reference_id].alignments[i]
                        == short_aln
                    ):
                        reference_maps_dict[short_aln.reference_id].alignments.pop(i)
                        break

            else:
                contig_map_label_start_containment, contig_map_label_end_containment = (
                    None,
                    None,
                )

                if label_1_start and label_1_end:
                    contig_map_label_start_containment = contig_maps_dict[
                        large_aln.map_id
                    ].get_label_position(
                        large_aln.get_corresponding_contig_map_label(label_1_start[0]),
                        large_aln.channel,
                    )
                    contig_map_label_end_containment = contig_maps_dict[
                        large_aln.map_id
                    ].get_label_position(
                        large_aln.get_corresponding_contig_map_label(label_1_end[0]),
                        large_aln.channel,
                    )
                elif label_2_start and label_2_end:
                    contig_map_label_start_containment = contig_maps_dict[
                        large_aln.map_id
                    ].get_label_position(
                        large_aln.get_corresponding_contig_map_label(label_2_start[0]),
                        large_aln.channel,
                    )
                    contig_map_label_end_containment = contig_maps_dict[
                        large_aln.map_id
                    ].get_label_position(
                        large_aln.get_corresponding_contig_map_label(label_2_end[0]),
                        large_aln.channel,
                    )
                else:
                    logging.debug(
                        f"Couldn't integrate map {short_aln.map_id} into {large_aln.map_id} (Anchor {large_aln.reference_id})because they don't share enough labels. Removing alignment of map {short_aln.map_id}"
                    )

                    for i in range(
                        0, len(reference_maps_dict[short_aln.reference_id].alignments)
                    ):
                        if (
                            reference_maps_dict[short_aln.reference_id].alignments[i]
                            == short_aln
                        ):
                            reference_maps_dict[short_aln.reference_id].alignments.pop(
                                i
                            )
                            break
                    return

                new_id = Key.get_max_id(key_dict) + 1

                logging.debug(
                    f"Map {large_aln.map_id} on contig {key_dict[(large_aln.map_id, large_aln.channel, large_aln.reference_id)][0]} will be broken in two maps."
                )
                logging.debug(
                    f"Map {large_aln.map_id} before change: {large_aln.reference_start} -> {large_aln.reference_end}"
                )

                new_map_start, new_map_end, new_map_size = 0, 0, 0
                old_map_start, old_map_end, old_map_size = 0, 0, 0
                new_aln_start, new_aln_end = 0, 0
                if large_aln.orientation == "+":
                    new_map_start = (
                        key_dict[
                            (
                                large_aln.map_id,
                                large_aln.channel,
                                large_aln.reference_id,
                            )
                        ][1]
                        + contig_map_label_end_containment
                    )
                    new_map_end = key_dict[
                        (large_aln.map_id, large_aln.channel, large_aln.reference_id)
                    ][2]
                    new_map_size = new_map_end - new_map_start
                    old_map_start = key_dict[
                        (large_aln.map_id, large_aln.channel, large_aln.reference_id)
                    ][1]
                    old_map_end = (
                        key_dict[
                            (
                                large_aln.map_id,
                                large_aln.channel,
                                large_aln.reference_id,
                            )
                        ][1]
                        + contig_map_label_start_containment
                    )
                    new_aln_start = 1
                    new_aln_end = large_aln.map_end - contig_map_label_end_containment
                else:
                    new_map_start = key_dict[
                        (large_aln.map_id, large_aln.channel, large_aln.reference_id)
                    ][1]
                    new_map_end = (
                        key_dict[
                            (
                                large_aln.map_id,
                                large_aln.channel,
                                large_aln.reference_id,
                            )
                        ][1]
                        + contig_map_label_end_containment
                    )
                    old_map_start = (
                        key_dict[
                            (
                                large_aln.map_id,
                                large_aln.channel,
                                large_aln.reference_id,
                            )
                        ][1]
                        + contig_map_label_start_containment
                    )
                    old_map_end = key_dict[
                        (large_aln.map_id, large_aln.channel, large_aln.reference_id)
                    ][2]
                    new_aln_end = 1
                    new_aln_start = contig_map_label_end_containment
                new_map_size = new_map_end - new_map_start
                old_map_size = old_map_end - old_map_start

                key_dict[(new_id, large_aln.channel, large_aln.reference_id)] = (
                    key_dict[
                        (large_aln.map_id, large_aln.channel, large_aln.reference_id)
                    ][0],
                    new_map_start,
                    new_map_end,
                    new_map_size,
                )
                key_dict[
                    (large_aln.map_id, large_aln.channel, large_aln.reference_id)
                ] = (
                    key_dict[
                        (large_aln.map_id, large_aln.channel, large_aln.reference_id)
                    ][0],
                    old_map_start,
                    old_map_end,
                    old_map_size,
                )

                contig_maps_dict[new_id] = Map.Map(
                    new_id,
                    contig_maps_dict[large_aln.map_id].labels_1,
                    contig_maps_dict[large_aln.map_id].labels_2,
                )

                new_aln = Alignment(
                    new_id,
                    large_aln.reference_id,
                    new_aln_start,
                    new_aln_end,
                    copy.deepcopy(short_aln.reference_end),
                    copy.deepcopy(large_aln.reference_end),
                    large_aln.channel,
                    large_aln.orientation,
                    large_aln.alignment,
                )
                reference_maps_dict[large_aln.reference_id].add_alignment(new_aln)
                new_aln.update_alignments(contig_maps_dict)

                if label_1_start and label_1_end:
                    large_aln.reference_end = label_1_start[1]
                elif label_2_start and label_2_end:
                    large_aln.reference_end = label_2_start[1]

                if large_aln.orientation == "+":
                    large_aln.map_end = contig_map_label_start_containment
                else:
                    large_aln.map_end = 1
                    large_aln.map_start -= contig_map_label_start_containment
                large_aln.update_alignments(contig_maps_dict)

                logging.debug(
                    f"Map {large_aln.map_id} after change: {large_aln.reference_start} -> {large_aln.reference_end}"
                )
                logging.debug(
                    f"Map {short_aln.map_id}: {short_aln.reference_start} -> {short_aln.reference_end}"
                )
                logging.debug(
                    f"Map {new_aln.map_id}: {new_aln.reference_start} -> {new_aln.reference_end}"
                )


def print_agp_line_with_intersection(
    aln_1,
    aln_2,
    previous_ref_end,
    previously_scaffolded_maps,
    contigs_map_dict,
    previous_part_number,
    key_dict,
):
    """
    Prints a line formatted by following the AGP2 standard. Used when two contig maps share labels.
    
    :param aln_1: Alignment that has the smallest reference_start
    :type aln_1: Alignment
    :param aln_2: Alignment that has the highest reference_start
    :type aln_2: Alignment
    :param previous_ref_end: Current position in the scaffold that is being built
    :type previous_ref_end: int
    :param previously_scaffolded_maps: Contig map ids that were previously used to build the current scaffold
    :type previously_scaffolded_maps: list(int)
    :param contigs_map_dict: Dict containing contig maps
    :type contigs_map_dict: dict(int, Map)
    :param previous_part_number: Id of the previous line
    :type previous_part_number: int
    :param key_dict: Dict containing the correspondance between contigs and contig maps
    :type key_dict: dict((int, int, int), (str, int, int, int))
    :return: Current position in the scaffold after applying the changes
    :rtype: int
    """

    intersection = get_shared_labels(aln_1, aln_2)
    label_1, label_2 = None, None
    label_1_pos, label_2_pos = None, None

    if intersection[0] or intersection[1]:
        try:
            label_1 = aln_1.get_corresponding_contig_map_label(
                intersection[aln_1.channel - 1][-1]
            )
            label_1_pos = contigs_map_dict[aln_1.map_id].get_label_position(
                label_1, aln_1.channel
            )
        except:
            logging.debug(f"Couldn't find label_1 on map {aln_1.map_id}")

        try:
            label_2 = aln_2.get_corresponding_contig_map_label(
                intersection[aln_2.channel - 1][-1]
            )
            label_2_pos = contigs_map_dict[aln_2.map_id].get_label_position(
                label_2, aln_2.channel
            )
        except:
            logging.debug(f"Couldn't find label_2 on map {aln_2.map_id}")

        if not label_1_pos or not label_2_pos:
            logging.debug(
                f"Couldn't detect the position of shared labels between maps {aln_1.map_id} and {aln_2.map_id} on anchor {aln_1.reference_id}"
            )
            intersection = None

        if intersection:
            start_map_1, start_map_2 = None, None
            end_map_1 = None

            # First contig map we encounter
            if aln_1.map_id not in previously_scaffolded_maps:
                if aln_1.orientation == "+":
                    start_map_1 = key_dict[
                        (aln_1.map_id, aln_1.channel, aln_1.reference_id)
                    ][1]
                else:
                    start_map_1 = key_dict[
                        (aln_1.map_id, aln_1.channel, aln_1.reference_id)
                    ][2]

            else:
                start_map_1 = previously_scaffolded_maps[aln_1.map_id][0]

            end_map_1 = (
                label_1_pos
                + key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][1]
            )
            start_map_2 = (
                label_2_pos
                + key_dict[(aln_2.map_id, aln_2.channel, aln_1.reference_id)][1]
            )

            if start_map_1 != end_map_1:
                with open("Phase_1/phase_1.agp", "a") as out:
                    out.write(
                        f"Super-Scaffold_{aln_1.reference_id}\t{int(previous_ref_end)}\t{int(previous_ref_end + math.fabs(end_map_1 - start_map_1))}\t{previous_part_number + 1}\tW\t{key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][0]}\t{min(start_map_1, end_map_1)}\t{max(end_map_1, start_map_1)}\t{aln_1.orientation}\n"
                    )

                previously_scaffolded_maps[aln_1.map_id] = (start_map_1, end_map_1)
                previously_scaffolded_maps[aln_2.map_id] = (start_map_2, None)
                previous_ref_end = (
                    previous_ref_end + math.fabs(end_map_1 - start_map_1) + 1
                )
                previous_part_number += 1
    else:
        intersection = None

    return (previous_ref_end, intersection, previous_part_number)


def print_gap_line(
    aln_1, aln_2, reference_id, previous_reference_end, previous_part_number, key_dict
):
    """
    Prints an AGP line, formatted following the AGP2 standard. Used to print an 'N' line.
    
    :param aln_1: Alignment that has the smallest reference_start
    :type aln_1: Alignment
    :param aln_2: Alignment that has the highest reference_start
    :type aln_2: Alignment
    :param reference_id: Id of the anchor map
    :type reference_id: int
    :param previous_reference_end: Current position in the scaffold being built
    :type previous_reference_end: int
    :param previous_part_number: Current AGP line id
    :type previous_part_number: int
    :param key_dict: Dict containing the correspondance between contigs and contig maps
    :type key_dict: dict((int, int, int), (str, int, int, int))
    :return: Current position in the scaffold being built after applying changes
    :rtype: int
    """

    gapsize = 0
    delta_1, delta_2, space_between_contigs = 0, 0, 0

    if aln_1.orientation == "+":
        delta_1 = (
            key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][3]
            - aln_1.map_end
        )
    else:
        delta_1 = aln_1.map_end

    if aln_2.orientation == "+":
        delta_2 = aln_2.map_start
    else:
        delta_2 = (
            key_dict[(aln_2.map_id, aln_2.channel, aln_2.reference_id)][3]
            - aln_2.map_start
        )

    space_between_contigs = (
        (aln_2.reference_start - aln_1.reference_end) - delta_1 - delta_2
    )

    if space_between_contigs <= 100:
        gapsize = 13
    else:
        gapsize = space_between_contigs

    with open("Phase_1/phase_1.agp", "a") as out:
        out.write(
            f"Super-Scaffold_{reference_id}\t{int(previous_reference_end)}\t{int(previous_reference_end + gapsize - 1)}\t{previous_part_number + 1}\tN\t{gapsize}\tscaffold\tyes\tmap\n"
        )

    return previous_reference_end + gapsize, previous_part_number + 1


def print_agp_line_no_intersection(
    aln_1,
    aln_2,
    previous_ref_end,
    previously_scaffolded_maps,
    key_dict,
    previous_part_number,
):
    """
    Prints an AGP line, formatted following the AGP2 standard. Used when two contig maps don't share anchor labels.
    
    :param aln_1: Alignment that has the smallest reference_start
    :type aln_1: Alignment
    :param aln_2: Alignment that has the highest reference_start
    :type aln_2: Alignment
    :param previous_ref_end: Current position in the scaffold being built
    :type previous_ref_end: int
    :param previously_scaffolded_maps: Contig map ids that were previously used to build the current scaffold
    :type previously_scaffolded_maps: list(int)
    :param key_dict: Dict containing the correspondance between contigs and contig maps
    :type key_dict: dict((int, int, int), (str, int, int, int))
    :param previous_part_number: Current AGP line id
    :type previous_part_number: int
    :return: Current position in the scaffold being build after applying the changes
    :rtype: int
    """

    start_map_1, start_map_2 = None, None
    end_map_1 = None

    if aln_1.map_id not in previously_scaffolded_maps:
        if aln_1.orientation == "+":
            start_map_1 = key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][1]
        else:
            start_map_1 = key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][2]
    else:
        start_map_1 = previously_scaffolded_maps[aln_1.map_id][0]

    if aln_1.orientation == "+":
        end_map_1 = key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][2]
    else:
        end_map_1 = key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][1]

    if aln_2:
        if aln_2.orientation == "+":
            start_map_2 = key_dict[(aln_2.map_id, aln_2.channel, aln_2.reference_id)][1]
        else:
            start_map_2 = key_dict[(aln_2.map_id, aln_2.channel, aln_2.reference_id)][2]
        previously_scaffolded_maps[aln_2.map_id] = (start_map_2, None)

    with open("Phase_1/phase_1.agp", "a") as out:
        out.write(
            f"Super-Scaffold_{aln_1.reference_id}\t{int(previous_ref_end)}\t{int(previous_ref_end + math.fabs(end_map_1 - start_map_1))}\t{previous_part_number + 1}\tW\t{key_dict[(aln_1.map_id, aln_1.channel, aln_1.reference_id)][0]}\t{min(start_map_1, end_map_1)}\t{max(end_map_1, start_map_1)}\t{aln_1.orientation}\n",
        )

    previous_ref_end = int(previous_ref_end + math.fabs(end_map_1 - start_map_1)) + 1
    previous_part_number += 1

    if aln_2:
        previous_ref_end, previous_part_number = print_gap_line(
            aln_1,
            aln_2,
            aln_1.reference_id,
            previous_ref_end,
            previous_part_number,
            key_dict,
        )

    previously_scaffolded_maps[aln_1.map_id] = (start_map_1, end_map_1)
    return previous_ref_end, previous_part_number


def print_agp(
    reference_maps_dict, key_dict, deleted_xmap_records, contigs_map_dict,
):
    """
    Searches for shared labels between two Alignment objects and calls the correct function 
    
    :param reference_maps_dict: Dict containing anchor maps
    :type reference_maps_dict: dict(int, Map)
    :param key_dict: Dict containing the correspondance between contigs and contig maps
    :type key_dict: dict((int, int, int): (str, int, int, int))
    :param deleted_xmap_records: Dict containing smaller alignments that weren't retained when parsing xmaps
    :type deleted_xmap_records: dict(integer: Alignment)
    :param contigs_map_dict: Dict containing contig maps
    :type contigs_map_dict: dict(int, Map)
    :return: List containing contig maps that were used to build current scaffold
    :rtype: list(int)
    """

    f = open("Phase_1/phase_1.agp", "w")
    f.close()
    scaffolded_maps = []
    for reference_map_id in reference_maps_dict:
        previous_ref_end = 1
        previously_scaffolded_maps = defaultdict(list)
        previous_part_number = 0
        intersection = None
        aln_1, aln_2 = None, None

        for i in range(0, len(reference_maps_dict[reference_map_id].alignments) - 1):
            intersection = None
            aln_1, aln_2 = None, None
            aln_1 = reference_maps_dict[reference_map_id].alignments[i]
            aln_2 = reference_maps_dict[reference_map_id].alignments[i + 1]

            same_channel_alignments_found = True

            if aln_1.channel != aln_2.channel:
                aln_1_change_gap = 1000000000000
                aln_2_change_gap = 1000000000000
                if aln_1.map_id in deleted_xmap_records:
                    aln_1_change_gap = (
                        aln_2.reference_start
                        - deleted_xmap_records[aln_1.map_id].reference_end
                    )
                if aln_2.map_id in deleted_xmap_records:
                    aln_2_change_gap = (
                        deleted_xmap_records[aln_2.map_id].reference_start
                        - aln_1.reference_end
                    )

                if (
                    aln_1.map_id not in deleted_xmap_records
                    and aln_2.map_id not in deleted_xmap_records
                ):
                    logging.debug(
                        f"Couldn't find alignment of map {aln_1.map_id} and {aln_2.map_id} on the same channel"
                    )
                    same_channel_alignments_found = False
                else:
                    if aln_1_change_gap <= aln_2_change_gap:
                        aln_1 = deleted_xmap_records[aln_1.map_id]
                    else:
                        aln_2 = deleted_xmap_records[aln_2.map_id]

            if same_channel_alignments_found and aln_1.reference_end > aln_2.reference_start:
                (
                    previous_ref_end,
                    intersection,
                    previous_part_number,
                ) = print_agp_line_with_intersection(
                    aln_1,
                    aln_2,
                    previous_ref_end,
                    previously_scaffolded_maps,
                    contigs_map_dict,
                    previous_part_number,
                    key_dict,
                )

            if not same_channel_alignments_found or not intersection:
                previous_ref_end, previous_part_number = print_agp_line_no_intersection(
                    aln_1,
                    aln_2,
                    previous_ref_end,
                    previously_scaffolded_maps,
                    key_dict,
                    previous_part_number,
                )

        if len(reference_maps_dict[reference_map_id].alignments) > 1:
            print_agp_line_no_intersection(
                aln_2,
                None,
                previous_ref_end,
                previously_scaffolded_maps,
                key_dict,
                previous_part_number,
            )
        elif len(reference_maps_dict[reference_map_id].alignments) == 1:
            print_agp_line_no_intersection(
                reference_maps_dict[reference_map_id].alignments[0],
                None,
                previous_ref_end,
                previously_scaffolded_maps,
                key_dict,
                previous_part_number,
            )

        scaffolded_maps.extend(previously_scaffolded_maps.keys())
    return scaffolded_maps


def write_unplaced_contigs(key_dict, contigs_sequence_dict, scaffolded_maps):
    """
    Incorporates contigs that weren't scaffolded into the AGP file
    
    :param key_dict: Dict containing the correspondance between contigs and contig maps
    :type key_dict: dict((int, int, int), (str, int, int, int))
    :param contigs_sequence_dict: Dict containing fasta sequences
    :type contigs_sequence_dict: dict(str, str)
    :param scaffolded_maps: List containing contig map ids that were scaffolded
    :type scaffolded_maps: list(int)
    """

    logging.info("Writing unplaced contigs in agp file")
    counter = 1
    with open("Phase_1/phase_1.agp", "a") as agp:
        for contig_map, channel in key_dict:
            if contig_map not in scaffolded_maps:
                contig_name = key_dict[(contig_map, channel)][0]
                contig_start = key_dict[(contig_map, channel)][1]
                contig_end = key_dict[(contig_map, channel)][2] - 1
                scaffold_name = f"{key_dict[(contig_map, channel)][0]}_subseq_{contig_start}:{contig_end}"
                scaffold_start = 1
                scaffold_end = contig_end - contig_start + 1
                orientation = "+"
                agp.write(
                    f"{scaffold_name}\t{scaffold_start}\t{scaffold_end}\t1\tW\t{contig_name}\t{contig_start}\t{contig_end}\t{orientation}\n"
                )
                counter += 1
                scaffolded_maps.append(contig_map)


def solve_alignment_containment(reference_maps_dict, contigs_map_dict, key_dict):
    """
    Calls the contained alignment solver function for each alignment couple
    
    :param contained_alignments: Tuple containing the contained alignment (second position) and the large alignment (first position)
    :type reference_maps_dict: dict(int, Map)
    :param contigs_map_dict: Dict containing contig maps
    :type contigs_map_dict: dict(int, Map)
    :param key_dict: Dict containing correspondance between contigs and contig maps
    :type key_dict: dict((int, int, int), (str, int, int, int))
    """

    for i in range(0, 4):
        Map.sort_map_alignments(reference_maps_dict)
        contained_alignments = Map.check_map_containment(reference_maps_dict)

        if contained_alignments:
            for i in range(0, len(contained_alignments)):
                for aln_couple in contained_alignments[i]:
                    solve_containment(
                        aln_couple, reference_maps_dict, contigs_map_dict, key_dict
                    )

        else:
            break
