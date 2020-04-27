from biscot import Misc

import copy
import logging


class Map:
    def __init__(self, map_id, labels_1, labels_2, alignments=[]):
        self.map_id = map_id
        self.labels_1 = labels_1
        self.labels_2 = labels_2
        self.nb_channel_1_labels = len(labels_1)
        self.nb_channel_2_labels = len(labels_2)
        self.alignments = copy.deepcopy(alignments)
        self.contig_maps = []

    def add_channel_1_label(self, label_position):
        """
        Adds a label to the list of channel 1 labels
        
        :param label_position: Position of the label on the map
        :type label_position: int
        """

        self.nb_channel_1_labels += 1
        self.labels_1.append((self.nb_channel_1_labels, label_position))

    def add_channel_2_label(self, label_position):
        """
        Adds a label to the list of channel 2 labels
        
        :param label_position: Position of the label on the map
        :type label_position: int
        """

        self.nb_channel_2_labels += 1
        self.labels_2.append((self.nb_channel_2_labels, label_position))

    def get_label_position(self, label_id, channel):
        """
        Returns a label position on a map based on its id and channel
        
        :param label_id: Label id to look for
        :type label_id: integer
        :param channel: Enzyme channel of the searched label
        :type channel: int
        :raises Exception: If the label couldn't be found
        :return: Searched label position
        :rtype: int
        """

        if channel == 1:
            for label, position in self.labels_1:
                if label == label_id:
                    return position
        elif channel == 2:
            for label, position in self.labels_2:
                if label == label_id:
                    return position
        raise Exception(
            f"Didn't find label {label_id} on map {self.map_id} on channel {channel}."
        )

    def add_alignment(self, aln):
        """
        Adds an Alignment object to the list of alignments
        
        :param aln: Alignment to add
        :type aln: Alignment
        """

        self.alignments.append(aln)
        self.contig_maps.append(aln.map_id)

    def sort_alignments(self):
        """
        Sorts the list of alignments
        """

        self.alignments = sorted(
            self.alignments, key=lambda alignment: alignment.reference_start
        )

    def check_containment(self):
        """
        Parses the list of alignments in search of alignments that could be contained into another one, i.e. reference_start_aln_1 < reference_start_aln_2 and reference_end_aln_1 > reference_end_aln_2
        
        :return: List of tuples containing the contained alignment at the second position and the alignment containing it at the first position
        :rtype: list(tuple(Alignment, Alignment))
        """

        contained_alns = []
        for i in range(0, len(self.alignments) - 1):
            aln_1 = self.alignments[i]

            for j in range(i + 1, len(self.alignments)):
                aln_2 = self.alignments[j]

                if (
                    aln_1.reference_start < aln_2.reference_start
                    and aln_1.reference_end > aln_2.reference_end
                ):
                    logging.debug(
                        f"Map {aln_2.map_id} ({aln_2.reference_start} -> {aln_2.reference_end}) is contained in map {aln_1.map_id} ({aln_1.reference_start} -> {aln_1.reference_end}) on anchor {aln_1.reference_id} channel {aln_1.channel}"
                    )
                    contained_alns.append((aln_1, aln_2))

                elif (
                    aln_2.reference_start < aln_1.reference_start
                    and aln_2.reference_end > aln_1.reference_end
                ):
                    logging.debug(
                        f"Map {aln_1.map_id} ({aln_1.reference_start} -> {aln_1.reference_end}) is contained in map {aln_2.map_id} ({aln_2.reference_start} -> {aln_2.reference_end}) on anchor {aln_2.reference_id} channel {aln_2.channel}"
                    )
                    contained_alns.append((aln_2, aln_1))

        return contained_alns

    def print_alignments(self):
        """
        Prints the alignments of a Map object
        """

        for aln in self.alignments:
            print(aln, flush=True)

    def __str__(self):
        txt = f"{self.map_id}\t{self.nb_channel_1_labels}\t{self.nb_channel_2_labels}\t{self.alignments}"
        return txt


def parse_reference_cmap(reference_cmap_file_path):
    """
    Parses a reference CMAP file to extract anchor labels
    
    :param reference_cmap_file_path: Path to a CMAP file
    :type reference_cmap_file_path: str
    :return: Dict containing anchor maps
    :rtype: dict(int, Map)
    """

    reference_cmap_file = open(reference_cmap_file_path)
    reference_maps_dict = {}

    for line in reference_cmap_file:
        if not line.startswith("#"):
            line = line.rstrip("\n").split("\t")

            map_id = int(line[0])
            label_channel = int(line[4])
            label_position = int(line[5].split(".")[0])

            if map_id not in reference_maps_dict:
                reference_maps_dict[map_id] = Map(map_id, [], [])

            if label_channel == 1:
                reference_maps_dict[map_id].add_channel_1_label(label_position)
            elif label_channel == 2:
                reference_maps_dict[map_id].add_channel_2_label(label_position)
    reference_cmap_file.close()

    return reference_maps_dict


def parse_contig_cmap(cmap_1_path, cmap_2_path):
    """
    Parses one or two contig CMAP files to extract contig labels
    
    :param cmap_1_path: Path to a CMAP file
    :type cmap_1_path: str
    :param cmap_2_path: Path to a CMAP file
    :type cmap_2_path: str
    :return: Dict containing contg maps
    :rtype: dict(str: Map)
    """

    cmap_1_file = open(cmap_1_path)
    contigs_map_dict = {}

    for line in cmap_1_file:
        if not line.startswith("#"):
            line = line.rstrip("\n").split("\t")

            map_id = int(line[0])
            label_position = int(line[5].split(".")[0])

            if map_id not in contigs_map_dict:
                contigs_map_dict[map_id] = Map(map_id, [], [])

            contigs_map_dict[map_id].add_channel_1_label(label_position)
    cmap_1_file.close()

    if cmap_2_path:
        cmap_2_file = open(cmap_2_path)

        for line in cmap_2_file:
            if not line.startswith("#"):
                line = line.rstrip("\n").split("\t")

                map_id = int(line[0])
                label_position = int(line[5].split(".")[0])

                if map_id not in contigs_map_dict:
                    contigs_map_dict[map_id] = Map(map_id, [], [])

                contigs_map_dict[map_id].add_channel_2_label(label_position)
        cmap_2_file.close()

    return contigs_map_dict


def sort_map_alignments(reference_maps_dict):
    """
    Sorts the alignments of a Map object by reference_start value
    
    :param reference_maps_dict: Dict containing anchor maps
    :type reference_maps_dict: dict(int, Map)
    """

    logging.info("Sorting alignments")
    for map in reference_maps_dict:
        reference_maps_dict[map].sort_alignments()


def check_map_containment(reference_maps_dict):
    """
    Parses all alignments of a map in search of contaned alignments
    
    :param reference_maps_dict: Dict containing anchor maps
    :type reference_maps_dict: dict(int: Map)
    :return: A list containing containing contained alignments
    :rtype: list(tuple(Alignment, Alignment))
    """

    logging.info("Looking for contained maps")
    contained_alignments = []
    for anchor in reference_maps_dict:
        alns = reference_maps_dict[anchor].check_containment()
        if len(alns) > 0:
            contained_alignments.append(alns)
    return contained_alignments
