import os
import argparse
import argcomplete
import sys
import gzip

def kmer_counting(seq: str, k: int) -> dict:
    """
    Counts the kmers in a sequence.
    :param seq: sequence to kmer count
    :param k: kmer length
    :return kmer_diz: dictionary with the kmers and their counts:
    """
    kmer_set = set([seq[i:i+k] for i in range(len(seq)-k+1)])
    kmer_diz = {kmer:seq.count(kmer) for kmer in kmer_set}
    return kmer_diz

def check_file(file_path: str) -> bool:
    try:
        stat_file = os.stat(file_path)
        res = True
    except FileNotFoundError:
        sys.exit(f"{file_path} does not exist!")
    else:
        return res

def findComp(sequence:str)->str:
    """
    function which returns complementary sequence of DNA
    input type: str
    param input: sequence you want the complementary one
    return type:str
    """
    cp = {"A":"T", "C":"G", "G":"C", "T":"A","R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K","B":"V","D":"H","H":"D","V":"B","N":"N"}
    comp = [cp[nt] for nt in sequence]
    comp = "".join(comp)
    return comp

def findCompRev(comp:str)->str:
    """
    function which returns complementary reverse sequence of DNA
    input type:str
    param type: sequence you want the complementary reverse one
    return type:str
    """
    comp_rev = comp[::-1]
    return comp_rev


def extract_info(filename,nuovo_file):
    if filename.endswith('.gz') or filename.endswith('.gzip'):
        with gzip.open(filename, 'rt') as file:
            contenuto = file.readlines()
            output_file = nuovo_file + '.fastq'
        with open(output_file, 'w') as output:
            for line in contenuto:
                output.write(line)

    elif filename.endswith('.fastq') or filename.endswith('.fq'):
        with open(filename, 'r') as file:
            c = file.readlines()
            output_file = nuovo_file + '.fastq'
        with open(output_file, 'w') as output:
            for line in c:
                output.write(line)
    else:
        print("Wrong File Format. Choose a FastQ or GZ file")

    return output_file

def hard_trimming(sequence:str, qualList: list[int], qs: int, len_s:int):
    """
    Ad un certo punto mettiamo la documentazione
    """
    # filter(lambda x,qs: x < qs, qualList)
    l = [x for x in qualList if x < qs]
    if l:
        i = qualList.index(l[0])
        qualList = qualList[:i]
        sequence = sequence[:i]
    if len(sequence) < len_s:
        sequence, qualList = None, None
    return sequence, qualList


def dinamic_trimming(quality: list[int], treshold=25, window=1):
    quality_inv = quality[::-1]
    list_mean = quality_inv[:window]
    if sum(list_mean)/len(list_mean) < treshold:
        for value in range(len(list_mean), len(quality_inv)):
            if (sum(list_mean)+quality_inv[value])/(len(list_mean)+1) < treshold:
                list_mean.append(quality_inv[value])
            else:
                break
        quality_inv_trim = quality_inv[len(list_mean):]
    else: quality_inv_trim = quality_inv
    return quality_inv_trim[::-1]