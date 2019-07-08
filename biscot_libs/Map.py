class Map :
    """
    A class used to represent a contig map

    Attributes :
        map_id : int
            A unique identifier
        contig_name : str
            A contig can be splitted in multiple parts. This attribute describes the name of the part that was used
        base_contig_name : str
            Name of the contig that was used for scaffolding BEFORE the split
        size : int
            Size of the map in nucleotides
        start : int
            Starting position of the map in the contig
        end :
            Ending position of the map in the contig
        labels : list((int, int, str))
            List of labels of the map. The tuple describes the label identifier, the label position and the label channel

    Methods :
        update_map(self, contig_name, dict_sequences)
            Updates the map coordinates and base contig name based on contig name
        add_label(self, cmap_splitted_line) 
            Appends a DLE or BspQI label to self.labels
        get_label_position_on_map(self, label)
            Returns the position of the specified label on the map
        get_label_position_on_contig(self, label, position)
            Returns the position of the specified label on the contig
    """

    def __init__(self, map_id, contig_name, size) :
        """
        Parameters :
            map_id : int
                Unique identifier
            contig_name : str
                Name of the contig after the split
            size : str
                Size of the map in nucleotides
        """

        self.map_id = map_id
        self.contig_name = contig_name
        self.base_contig_name = ""
        self.size = int(size)
        self.start = 0
        self.end = 0
        self.labels = []


    def update_map(self, contig_name, dict_sequences) :
        """
        Updates the map coordinates and base contig name based on contig name

        Parameters :
            contig_name : str
                Name of the contig the map originates from
            dict_sequences : dict{str: str}
                Dictionnary containing the name of the contigs as keys and their fasta sequences as values
        """

        start = -1
        end = -1

        # If the contig was splitted by the Bionano Access, we need to recompute coordinates
        if "subseq" in contig_name :
            contig_name_split = contig_name.split("_subseq_")
            self.start = min(int(contig_name_split[1].split(":")[0]), int(contig_name_split[1].split(":")[1])) - 1
            self.end = max(int(contig_name_split[1].split(":")[0]), int(contig_name_split[1].split(":")[1]))
            self.base_contig_name = contig_name_split[0]
            self.contig_name = contig_name

        else :
            self.start = 0
            self.end = len(dict_sequences[self.contig_name])
            self.base_contig_name = contig_name
            self.contig_name = contig_name


    def add_label(self, cmap_splitted_line, channel) :
        """
        Appends a DLE or BspQI label to self.labels

        Parameters :
            cmap_splitted_line : list(str)
                List of all fields of one CMAP line
        """

        label_id, label_position, label_channel = int(cmap_splitted_line[3]), int(cmap_splitted_line[5].split(".")[0]), int(cmap_splitted_line[4])
        if label_channel != 0 :
            self.labels.append((label_id, label_position, channel))


    def update_labels(self, difference) :
        """
        Removes labels that have a position > self.end or < self.start

        Parameters :
            difference : int
                Number to subtract from the position
        """

        tmp_labels = []
        for label_id, position, channel in self.labels :
            tmp_labels.append((label_id, position - difference, channel))
        self.labels = tmp_labels


    def get_label_position_on_map(self, label) :
        """
        Returns the position of the specified label on the map

        Parameters :
            label : int
                Identifier of the label to search for

        Returns :
            position : int
                Position of the specified label
        """

        for lab in self.labels :
            if lab[0] == label :
                return lab[1]


    def get_label_position_on_contig(self, label, position) :
        """
        Returns the position of the specified label on the contig

        Parameters :
            label : int 
                Identifier of the label to search for
            position : int
                Last known mapping position of the contig map

        Returns :
            position : list(int, int)
                Position of the DLE and BspQ1 labels on the contig map
        """
        positions = []

        # Internally, the Bionano access renames the DLE labels independently of the BspQI labels
        # We need to first parse DLE labels and then BspQI labels to account for this

        nb_label = 1
        for lab in self.labels :
            if lab[2] == 1 :
                if nb_label == label :
                    positions.append(self.start + lab[1])
                    break
                nb_label += 1

        nb_label = 1
        for lab in self.labels :
            if lab[2] == 2 :
                if nb_label == label :
                    positions.append(self.start + lab[1])
                    break
                nb_label += 1

        # In some cases, the label is not indicated in the CMAP
        # If so, returns the end mapping position of the contig map as it is always at a label
        if not positions :
            return [self.start + position]
        else :
            return positions


    def __str__(self) :
        txt = "%s\t%s\t%s\t%s\t%s\t%s" % (self.map_id, self.contig_name, self.base_contig_name, self.start, self.end, self.size)
        return txt
