import os
import argparse
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

def findCompRev(sequence:str)->str:
    """
    function which returns complementary sequence of DNA
    input type: str
    param input: sequence you want the complementary one
    return type:str
    """
    cp = {"A":"T", "C":"G", "G":"C", "T":"A","R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K","B":"V","D":"H","H":"D","V":"B","N":"N"}
    comp = [cp[nt] for nt in sequence]
    comp = "".join(comp)
    comp_rev = comp[::-1]
    return comp,comp_rev

def fn_comp_rev (nameseq):
    nameseq = nameseq.replace('\n','')
    reverse= nameseq[::-1]
    table= str.maketrans('ACTG','TGAC')
    comp = nameseq.translate(table)
    comp_rev= nameseq.translate(table)[::-1]
    return comp,comp_rev




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

