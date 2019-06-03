from biscot_libs import Anchor
import logging

class Alignment :
    """
    Class used to represent an alignment between a map and an anchor

    Atributes :
        line : str
            One line of the XMAP file that was used to create the Alignment object
        map_id : int 
            Identifier of the map
        map_start : int
            Start of the alignment in the map
        map_end : int
            End of the alignment in the map
        anchor_id : int
            Identifier of the anchor
        anchor_start : int
            Start of the alignment in the anchor
        anchor_end : int
            End of the alignment in the anchor
        orientation : str
            Strand of the alignment ('+' or '-')
        label_channel : str
            Defines which enzyme was used to obtain this label (typically '1' or '2')
        label_mappings : dict(int : list(int))
            Contains all map labels (values) that mapped against an anchor label (keys)

    Methods :
        get_anchor_labels(self) 
            Returns all anchor labels of an alignment
        get_corresponding_contig_map_label(self, anchor_label)
            Returns the first map label that was aligned against an anchor label in an alignment
        set_anchor_start(self, start)
            Modifies the anchor_start attribute of an alignment
        set_anchor_end(self, end)
            Modifies the anchor_end attribute of an alignment
        set_map_start(self, start)
            Modifies the map_start attribute of an alignment
        set_map_end(self, end)
            Modifies the map_end attribute of an alignment
        update_mappings(self, anchor)
            Changes the self.label_mappings attribute to reflect changes done to its anchor_start and anchor_end by removing some label alignments that are out of bounds
        find_shared_labels(anchor, map_id_1, map_id_2)
            Returns anchor labels that are shared in two different alignments (map overlap)
        get_genomic_intersection(mapping_pos_1, mapping_pos_2)
            Returns the overlap size between two alignment coordinates
    """

    def __init__(self, xmap_line) :
        """
        Parameters :
            xmap_line : str
                One line of a xmap file 
        """

        self.line = xmap_line
        xmap_line = xmap_line.strip().split("\t")

        self.map_id, self.map_start, self.map_end = int(xmap_line[1]), int(xmap_line[3].split(".")[0]), int(xmap_line[4].split(".")[0])
        self.anchor_id, self.anchor_start, self.anchor_end = int(xmap_line[2]), int(xmap_line[5].split(".")[0]), int(xmap_line[6].split(".")[0])
        self.orientation = xmap_line[7]
        self.label_channel = xmap_line[12]

        # Base format : "(anchor_label_1,map_label_1)(anchor_label_2, map_label_2)..."
        # New format  : ["anchor_label_1,map_label_1", "anchor_label_2,map_label_2"]
        label_mappings = xmap_line[13].replace("(", "").replace(")", " ").split(" ")

        self.label_mappings = {}
        for mapping in label_mappings :
            if mapping :
                anchor_label, map_label = mapping.split(",")
                try :
                    self.label_mappings[int(anchor_label)].append(int(map_label))
                except :
                    self.label_mappings[int(anchor_label)] = [int(map_label)]


    def get_anchor_labels(self) :
        """
        Returns all anchor labels of an alignment
        
        Returns :
            labels : set(int)
                List of anchor labels contained in the self.label_mappings attribute
        """

        return set(self.label_mappings.keys())


    def get_corresponding_contig_map_label(self, anchor_label) :
        """
        Returns the first map label that was aligned against an anchor label

        Parameters :
            anchor_label : int
                Anchor label to look for in the self.label_mappings attribute

        Returns :
            map_label : int
                First map label that was aligned against the specified anchor label
        """

        return self.label_mappings[anchor_label][0]


    def set_anchor_start(self, start) :
        """
        Modifies the anchor_start attribute of an alignment

        Parameters :
            start : int
                self.anchor_start will be changed to start
        """

        self.anchor_start = start


    def set_anchor_end(self, end) :
        """
        Modifies the anchor_end attribute of an alignment

        Parameters :
            end : int
                self.anchor_end will be changed to end
        """

        self.anchor_end = end


    def set_map_id(self, map_id) :
        """
        Modifies the map_map_id attribute of an alignment

        Parameters :
            start : int
                self.map_id will be changed to map_id
        """

        self.map_id = map_id


    def set_map_start(self, start) :
        """
        Modifies the map_start attribute of an alignment

        Parameters :
            start : int
                self.map_start will be changed to start
        """

        self.map_start = start


    def set_map_end(self, end) :
        """
        Modifies the map_end attribute of an alignment

        Parameters :
            end : int
                self.map_end will be changed to end
        """

        self.map_end = end


    def update_mappings(self, anchor) :
        """
        Changes the self.label_mappings attribute to reflect changes done to its anchor_start and anchor_end by removing some label alignments that are out of bounds

        Parameters :
            anchor : Anchor
                Anchor implied in the alignment

        Returns
            labels : tuple(int, int)
                First and last labels that have been removed
        """

        labels_to_remove = []

        # Don't remove the first n labels in order to still be able to chain maps after removing alignments
        nb_labels_to_keep = 0

        for i, anchor_label in enumerate(self.label_mappings) :
            try :
                if anchor.labels_DLE[anchor_label] > self.anchor_end or anchor.labels_DLE[anchor_label] < self.anchor_start :
                    if nb_labels_to_keep > 0 :
                        nb_labels_to_keep -= 1
                    else :
                        labels_to_remove.append(anchor_label)
            except :
                if anchor.labels_BspQI[anchor_label] > self.anchor_end or anchor.labels_BspQI[anchor_label] < self.anchor_start :
                    if nb_labels_to_keep > 0 :
                        nb_labels_to_keep -= 1
                    else :
                        labels_to_remove.append(anchor_label)

        for label in labels_to_remove :
            del self.label_mappings[label]

        logging.debug("Removed labels from %s to %s in alignment of map %s on anchor %s" % (labels_to_remove[0], labels_to_remove[-1], self.map_id, self.anchor_id))

        return (labels_to_remove[0], labels_to_remove[-1])


    def __str__(self) :
        txt = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.anchor_id, self.anchor_start, self.anchor_end, self.map_id, self.map_start, self.map_end, self.orientation)
        return txt


def find_shared_labels(anchor, map_id_1, map_id_2) :
    """
    Returns anchor labels that are shared in two different alignments (map overlap)

    Parameters :
        anchor : Anchor
            Anchor in which to search alignments
        map_id_1 : int 
            First map identifier to look for
        map_id_2 : int
            Second map identifier to look for

    Returns :
        intersection : set(int)
            Anchor labels that are shared by two alignments
    """

    mapping_dict = {}
    for aln in anchor :
        mapping_dict[aln.map_id] = aln.get_anchor_labels()

    intersection = mapping_dict[map_id_1].intersection(mapping_dict[map_id_2])
    intersection = sorted(intersection)
    return intersection


def get_overlap_size(mapping_pos_1, mapping_pos_2) :
    """
    Returns the overlap size between two alignment coordinates

    Parameters :
        mapping_pos_1 : tuple(int, int)
            Contains the coordinates of the first alignment
        mapping_pos_2 : tuple(int, int)
            Contains the coordinates of the second alignment

    Returns :
        intersection : int
            Size of the overlap between two alignments
    """

    start_1 = min(mapping_pos_1[0], mapping_pos_1[1])
    end_1 = max(mapping_pos_1[0], mapping_pos_1[1])
    start_2 = min(mapping_pos_2[0], mapping_pos_2[1])
    end_2 = max(mapping_pos_2[0], mapping_pos_2[1])

    return max(0, min(end_1, end_2) - max(start_1, start_2))


def parse_blat(scaffold_name, start) :
    """
    Parses the blat output file

    Parameters :
        scaffold_name : str
            Name of the scaffold where the gap is localized
        start : int
            Position of the gap
    
    Returns :
        reference_size : int
            Size of the query sequence
        aln_end : int
            Start of the alignment in reference
        aln_end : int
            End of the alignment in query
    """

    with open("blat_output.tmp") as blat :
        for i in range(0, 5) : blat.readline()

        max_score = 0
        max_score_line = ""

        for line in blat :
            line = line.split("\t")

            if int(line[0]) > max_score :
                max_score = int(line[0])
                max_score_line = line

        try :
            logging.debug("MATCH : %s %s %s", scaffold_name, start, max_score_line[0])
            return (int(max_score_line[14]), int(max_score_line[16]), int(max_score_line[12]))
        except :
            logging.debug("NO MATCH : %s %s" % (scaffold_name, start))
            return (None, None, None)
