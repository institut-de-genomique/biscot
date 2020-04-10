import copy
import logging


def parse_key(key_file_path):
    """
    Parses a Bionano '.key' file and extracts informations about contigs and their contig map counterparts
    
    :param key_file_path: Path to a '.key' file
    :type key_file_path: str
    :return: Dict containing the correspondance between contigs and contig maps
    :rtype: dict((id: int, channel: int), (contig_name: str, start: int, end:int, size: int))
    """

    key_dict = {}
    key_file = open(key_file_path)
    contig_name, start, end, size = None, None, None, None

    # Remove headers
    for line in key_file:
        if "CompntId" in line:
            break

    for line in key_file:
        component_id, component_name, component_length = line.rstrip("\n").split("\t")

        if "subseq" in component_name:
            contig_name, start_and_end = component_name.split("_subseq_")
            start, end = start_and_end.split(":")
        else:
            contig_name = component_name
            start = 1
            end = component_length

        component_id, start, end = int(component_id), int(start), int(end)
        size = end - start + 1

        # Create one key entry per channel
        key_dict[(component_id, 1)] = (contig_name, start, end, size)
        key_dict[(component_id, 2)] = (contig_name, start, end, size)

        if size != int(component_length):
            logging.info(
                f"WARNING: Map {component_id} (contig {contig_name}) has a wrong size"
            )

    return key_dict


def extend_key_dict(key_dict, reference_maps_dict):
    """
    Adds the reference id to the key_dict key as a contig can be placed multiple times and we don't want to modify its key informations erroneously in case of contained alignments
    
    :param key_dict: Dict containing the correspondance between contigs and contig maps
    :type key_dict: dict((int, int), (str, int, int, int))
    :param reference_maps_dict: Dict containing reference anchor maps
    :type reference_maps_dict: dict(int, Map)
    :return: Key dict containing the correspondance between contigs and contig maps, with the added information of reference id
    :rtype: dict((int, int, int), (str, int, int, int))
    """

    extended_key_dict = {}
    for reference_map in reference_maps_dict:
        for alignment in reference_maps_dict[reference_map].alignments:
            extended_key_dict[(alignment.map_id, 1, reference_map)] = copy.deepcopy(
                key_dict[(alignment.map_id, 1)]
            )
            extended_key_dict[(alignment.map_id, 2, reference_map)] = copy.deepcopy(
                key_dict[(alignment.map_id, 2)]
            )
    return extended_key_dict


def get_max_id(key_dict):
    """
    Gets the max id found inside a key_dict keys
    
    :param key_dict: Dict containing the correspondance between contigs and contig maps
    :type key_dict: dict((int, int, int), (str, int, int, int))
    :return: Maximum value of the key_dict keys
    :rtype: int
    """

    max_id = 0
    for map_id, _, _ in key_dict:
        if max_id < map_id:
            max_id = map_id
    return max_id
