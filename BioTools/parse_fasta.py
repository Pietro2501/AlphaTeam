import os
import sys
sys.path.append("C:/Users/Melania/Desktop/program_vsc/biotools")
#os.chdir('C:/Users/Melania/Desktop/program_vsc/biotools')
# biotools.tools import check_file
from Esame.tools1 import check_file

def parse_fasta(fileFasta: str) -> dict:
    """
    Parse a fasta file and return a dictionary containing che association among accession number and biological sequence.
    :param fileFasta: absolute or relative path to fasta file
    :type fileFasta: str
    :return seq_diz: dictionary if the file exists or otherwise None
    """
    seq_diz = None
    if check_file(fileFasta):
        #print(f"Il file {fileFasta} esiste e lo possiamo analizzare!!!")
        seq_diz = {}
        with open(fileFasta) as fasta:
            acc, seq = None, None
            line = fasta.readline().strip()
            while line:
                if line.startswith(">") and acc is None:
                    acc = line.split()[0][1:]
                    seq = ""
                elif line.startswith(">") and acc is not None:
                    seq_diz[acc] = seq
                    acc = line.split()[0][1:]
                    seq = ""
                else:
                    seq += line
                line = fasta.readline().strip()
            else:
                seq_diz[acc] = seq
    else:
        print(f"Il file {fileFasta} non esiste!!!")
    return seq_diz