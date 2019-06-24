import os
import logging
import sys
import subprocess
from collections import defaultdict
from collections import OrderedDict


_orig_print = print
def print(*args, **kwargs):
    """Unbuffured print function"""
    
    _orig_print(*args, flush=True, **kwargs)


def check_path(file_path) :
    """
    Checks if the file given as an argument exists

    Parameters :
        file_path : str
            File path to check

    Returns :
        bool
            Returns True if the file exists

    Raises :
        FileNotFoundError
            If the path doesn't exist
    """

    if not os.path.exists(file_path) :
        raise FileNotFoundError("%s file not found.")


def setup_logging(debug) :
    """
    Setups the logging streams

    Parameters :
        debug : bool
            True to activate debug logs
    """

    if debug :
        logging.basicConfig(format="[%(asctime)s]\t[%(filename)9.9s - %(lineno)03d]\t[%(levelname)5.5s]\t%(message)s", level=logging.DEBUG, filename="biscot.log")
    else :
        logging.basicConfig(format="[%(asctime)s]\t[%(filename)9.9s - %(lineno)03d]\t[%(levelname)5.5s]\t%(message)s", level=logging.INFO, filename="biscot.log")

    log_formatter = logging.Formatter("[%(asctime)s]\t[%(filename)9.9s - %(lineno)03d]\t[%(levelname)5.5s]\t%(message)s")
    root_logger = logging.getLogger()

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)


def launch_blat() :
    """Launches the blat alignment command during phase 2"""

    blat_cmd = "blat seq_1.tmp seq_2.tmp blat_output.tmp -minScore=3000"
    popen = subprocess.Popen(blat_cmd, shell = True, stdin = None, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = popen.communicate()
    #if err :
    #    logging.debug(err)


def mute_agp_file(agp_file_path, changes_file_path, prefix) :
    """
    Applies changes to an AGP file based on a file listing changes to be applied

    Parameters :
        agp_file_path : str
            Path to an AGP file
        changes_file_path : str
            Path to a file listing changes to be applied
        prefix : str
            Output agp prefix name
    """

    agp_file = open(agp_file_path)
    changes_file = open(changes_file_path)

    changes_dict = defaultdict(list)
    for line in changes_file :
        line = line.strip().split("\t")
        changes_dict[line[0]].append((line[1], line[2], line[3], line[4]))
    changes_file.close()

    agp_dict = OrderedDict()
    for line in agp_file :
        line = line.split("\t")
        try :
            agp_dict[line[0]].append(line)
        except :
            agp_dict[line[0]] = []
            agp_dict[line[0]].append(line)

    new_agp_lines = []
    # Mutating lines with a change 
    for scaffold in agp_dict :
        for i, line in enumerate(agp_dict[scaffold]) :
            if scaffold in changes_dict :
                change_found = False
                for change in changes_dict[scaffold] :
                    if change[0] == line[1] :
                        # Change start of next line on scaffold
                        try :
                            agp_dict[scaffold][i+1][1] = str(int(agp_dict[scaffold][i+1][1]) + int(change[2]))
                        except :
                            logging.info(change)
                            logging.info(agp_dict[scaffold][i+1])
                            exit()

                        # Change start of next line on contig
                        try :
                            agp_dict[scaffold][i+1][6] = str(int(agp_dict[scaffold][i+1][6]) + int(change[2]))
                        except:
                            logging.info(change)
                            logging.info(agp_dict[scaffold][i+1])
                            exit()

                        # Change end of previous line on scaffold
                        new_agp_lines[-1][2] = str(int(agp_dict[scaffold][i-1][2]) - (int(change[3]) - int(change[1])))

                        # Change end of previous line on contig
                        new_agp_lines[-1][7] = str(int(agp_dict[scaffold][i-1][7]) - (int(change[3]) - int(change[1])))

                        change_found = True
                        logging.debug("Removed gap at position %s on scaffold %s" % (line[1], scaffold))
                        break

                if not change_found :
                    new_agp_lines.append(agp_dict[scaffold][i])

            else :
                new_agp_lines.append(agp_dict[scaffold][i])

    with open(prefix + ".agp", "w") as out :
        last_scaffold = ""
        last_position = 0
        last_id = 0

        for line in new_agp_lines :
            if line[0] != last_scaffold :
                last_scaffold = ""
                last_position = 0
                last_id = 0

            if int(line[1]) != last_position + 1 :
                # New end = end - (start - last_pos)
                line[2] = str(int(line[2]) - (int(line[1]) - last_position))

                # New start = last_pos + 1
                line[1] = str(last_position + 1)

            line[3] = str(int(last_id) + 1)

            last_position = int(line[2])
            last_scaffold = line[0]
            last_id = line[3]

            out.write("\t".join(line))
