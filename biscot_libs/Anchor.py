from biscot_libs import Alignment
from collections import OrderedDict

class Anchor :
    """
    A class used to represent an hybrid scaffold used as an anchor by the Bionano Access

    Attributes :
        anchor_id : int 
            Identifier of the anchor
        alignments : list(Alignment)
            List of all alignments between maps and an anchor
        labels_DLE : OrderedDict
            Contains an anchor DLE labels as keys and their respective position as values
        labels_BspQI : OrderedDict
            Contains an anchor BspQI labels as keys and their respective position as values
        maps : list(int)
            Contains identifiers of maps that were aligned to an anchor

    Methods :
        add_alignment(self, alignment)
            Returns informations about the alignment between a map and an anchor
        sort_alignments(self, )
            Sorts the self.alignments attribute based on the anchor_start attribute of each alignment
        get_alignment_positions(self, contig_map_id)
            Returns alignment coordinates of the map 'contig_map_id' on the anchor
        add_DLE_label(self, cmap_splitted_line)
            Adds one DLE label entry in self.labels_DLE
        add_BspQI_label(self, cmap_splitted_line)
            Adds one BspQI label entry in self.labels_BspQI
        find_label_on_contig_map(self, contig_map_id, label)
            Returns the contig map label number of the map 'contig_map_id' that was aligned to the specified anchor label
        print_labels(self)
            Prints the label_id and corresponding position of the anchor DLE and BspQI labels
    """

    def __init__(self, anchor_id) :
        """
        Parameters :
            anchor_id : int
                A unique identifier corresponding to the super_scaffold number generated by the Bionano Access
        """

        self.anchor_id = anchor_id
        self.alignments = []
        self.labels_DLE = OrderedDict()
        self.labels_BspQI = OrderedDict()
        self.maps = []


    def add_alignment(self, alignment) :
        """
        Appends an Alignment object to the self.alignments attribute

        Parameters :
            alignment : Alignment
                Alignment object describing an alignment between a map and an anchor
        """

        self.alignments.append(alignment)


    def sort_alignments(self) :
        """
        Sorts the self.alignments attribute based on the anchor_start attribute of each Alignment object
        """

        self.alignments = sorted(self.alignments, key=lambda alignment: alignment.anchor_start)
        for aln in self.alignments :
            self.maps.append(aln.map_id)


    def get_alignment_positions(self, contig_map_id) :
        """
        Returns informations about the alignment between a map and an anchor

        Parameters :
            contig_map_id : int
                Unique identifier of a map 

        Returns :
            tuple(int, int, str, int, int) 
                Start of the alignment on anchor
                End of the alignment on anchor
                Orientation or strand of the alignment ('+' or '-')
                Start of the alignment on map
                End of the alignment on map
        """

        for aln in self.alignments :
            if aln.map_id == contig_map_id :
                return (aln.anchor_start, aln.anchor_end, aln.orientation, aln.map_start, aln.map_end)


    def add_DLE_label(self, cmap_splitted_line) :
        """
        Adds one DLE label entry in self.labels_DLE

        Parameters :
            cmap_splitted_line : list(string)
                list containing each field of one line of a CMAP file
        """

        label_id, label_position = int(cmap_splitted_line[3]), int(cmap_splitted_line[5].split(".")[0])
        self.labels_DLE[label_id] = label_position


    def add_BspQI_label(self, cmap_splitted_line) :
        """
        Adds one BspQI label entry in self.labels_BspQI

        Parameters :
            cmap_splitted_line : list(string)
                List containing each field of one line of a CMAP file
        """

        label_id, label_position = int(cmap_splitted_line[3]), int(cmap_splitted_line[5].split(".")[0])
        self.labels_BspQI[label_id] = label_position


    def find_label_on_contig_map(self, contig_map_id, label) :
        """
        Returns the contig map label number of the map 'contig_map_id' that was aligned to the specified anchor label

        Parameters :
            contig_map_id : int
                Map identifier
            label : int
                Label identifier on the anchor

        Returns :
            aln.get_corresponding_contig_map_label(label) : int
        """

        for aln in self.alignments :
            if aln.map_id == contig_map_id :
                try :
                    return aln.get_corresponding_contig_map_label(label)
                except :
                    continue


    def print_labels(self) :
        """Prints the label_id and corresponding position of the anchor DLE and BspQI labels"""

        for label in self.labels_DLE :
            print(label, self.labels_DLE[label])
        for label in self.labels_BspQI :
            print(label, self.labels_BspQI[label])


    def __iter__(self) :
        return iter(self.alignments)


    def __next__(self) :
        for aln in self.alignments :
            yield aln
