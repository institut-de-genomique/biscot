from biscot_libs import Anchor
import copy
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
        
        if xmap_line == "" :
            self.map_id = ""
            self.map_start = ""
            self.map_end = ""
            self.anchor_id = ""
            self.anchor_start = ""
            self.anchor_end = ""
            self.orientation = ""
            self.label_channel = ""
            self.label_mappings = ""

        else :
            xmap_line = xmap_line.strip().split("\t")

            self.map_id, self.map_start, self.map_end = int(xmap_line[1]), int(xmap_line[3].split(".")[0]), int(xmap_line[4].split(".")[0])
            self.anchor_id, self.anchor_start, self.anchor_end = int(xmap_line[2]), int(xmap_line[5].split(".")[0]), int(xmap_line[6].split(".")[0])
            self.orientation = xmap_line[7]
            self.label_channel = xmap_line[12]

            if self.anchor_start > self.anchor_end :
                self.anchor_start, self.anchor_end = self.anchor_end, self.anchor_start

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


    def get_anchor_labels_in_interval(self, start, end, anchor) :
        """
        Returns all anchor labels that were aligned to a map and that are between two positions

        Parameters :
            start : int
                Lower bound of positions
            end : int
                Upper bound of positions
            anchor : Anchor
                Anchor in which to search labels
 
        Retuns :
            labels_in_interval :
                List of mapped anchor labels between the upper bound and lower bound
        """

        labels_in_interval = []
        for label in self.label_mappings :
            try :
                if start < anchor.labels_DLE[label] < end :
                    labels_in_interval.append(label)
            except :
                try :
                    if start < anchor.labels_BspQI[label] < end :
                        labels_in_interval.append(label)  
                except :
                    pass
                pass

        return labels_in_interval


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


    def is_contained(self, map_1, aln_2, map_2) :
        """
        Verifies if map_1 is contained into map_2

        Parameters :
            map_1 : Map
            aln_2 : Alignment
                Alignment object of the second map
            map_2 : Map

        Returns :
            bool
                True if contained, False otherwise
        """

        aln_1 = self
        start_containment = aln_1.anchor_start >= aln_2.anchor_start
        end_containment = aln_1.anchor_end <= aln_2.anchor_end

        if map_2.map_id == 3697 and map_1.map_id == 81:
            print(map_1.map_id, map_2.map_id)
            print(start_containment)
            print(end_containment)

#        if aln_1.orientation == "+" and aln_2.orientation == "+" :
#            start_containment = aln_1.anchor_start - aln_1.map_start > aln_2.anchor_start - aln_2.map_start
#            end_containment = aln_1.anchor_end + (map_1.size - aln_1.map_end) < aln_2.anchor_end + (map_2.size - aln_2.map_end)
#
#        elif aln_1.orientation == "-" and aln_2.orientation == "-" :
#            start_containment = aln_1.anchor_start - aln_1.map_end > aln_2.anchor_start - aln_2.map_end
#            end_containment = aln_1.anchor_end + (map_1.size - aln_1.map_start) < aln_2.anchor_end + (map_2.size - aln_2.map_start)
#
#        elif aln_1.orientation == "+" and aln_2.orientation == "-" :
#            start_containment = aln_1.anchor_start - aln_1.map_start > aln_2.anchor_start - aln_2.map_end
#            end_containment = aln_1.anchor_end + (map_1.size - aln_1.map_end) < aln_2.anchor_end + (map_2.size - aln_2.map_start)
#
#        elif aln_1.orientation == "-" and aln_2.orientation == "+" :
#            start_containment = aln_1.anchor_start - aln_1.map_end > aln_2.anchor_start - aln_2.map_start
#            end_containment = aln_1.anchor_end + (map_1.size - aln_1.map_start) < aln_2.anchor_end + (map_2.size - aln_2.map_end)

        return(start_containment and end_containment)


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
            labels : tuple(tuple(int, int), tuple(int, int))
                First and last labels that have been removed and the corresponding labels on map
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
                        labels_to_remove.append((anchor_label, self.label_mappings[anchor_label][0]))
            except :
                if anchor.labels_BspQI[anchor_label] > self.anchor_end or anchor.labels_BspQI[anchor_label] < self.anchor_start :
                    if nb_labels_to_keep > 0 :
                        nb_labels_to_keep -= 1
                    else :
                        labels_to_remove.append((anchor_label, self.label_mappings[anchor_label][0]))

        for label in labels_to_remove :
            del self.label_mappings[label[0]]

        try :
            logging.debug("Removed labels from %s to %s in alignment of map %s on anchor %s" % (labels_to_remove[0], labels_to_remove[-1], self.map_id, self.anchor_id))
        except :
            logging.debug("Removed 0 label, no mappings at this position?")

        try :
            return (labels_to_remove[0], labels_to_remove[-1])
        except :
            return (None, None)


    def __str__(self) :
        txt = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.anchor_id, self.anchor_start, self.anchor_end, self.map_id, self.map_start, self.map_end, self.orientation)
        return txt


def copy_alignment(aln_2) :
    """
    Copy all paramters from a map alignment into a new Alignment object

    Parameters :
        aln_2 : Alignment
            Alignment object to copy parameters from
    """

    tmp = Alignment("")

    tmp.map_id = aln_2.map_id
    tmp.map_start = aln_2.map_start
    tmp.map_end = aln_2.map_end
    tmp.anchor_id = aln_2.anchor_id
    tmp.anchor_start = aln_2.anchor_start
    tmp.anchor_end = aln_2.anchor_end
    tmp.orientation = aln_2.orientation
    tmp.label_channel = aln_2.label_channel
    tmp.label_mappings = copy.deepcopy(aln_2.label_mappings)

    return(tmp)


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
