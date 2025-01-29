import os
import sys
import gzip
import pandas as pd

def check_file(file_path: str) -> bool:
    try:
        stat_file = os.stat(file_path)
        res = True
    except FileNotFoundError:
        sys.exit(f"{file_path} does not exist!")
    else:
        return res

def fn_comp_rev(nameseq):
    nameseq = nameseq.replace('\n','')
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

def hard_trimming(sequence:str, qualList: list[int], qs: int, len_s:int):
    l = [x for x in qualList if x < qs]
    if l:
        i = qualList.index(l[0])
        qualList = qualList[:i]
        sequence = sequence[:i]
    if len(sequence) < len_s:
        sequence, qualList = None, None
    return sequence, qualList

def dynamic_trimming(seq: str, qscore: list[int], threshold_q: int = 20, min_length: float = 0.95) -> tuple[str, list[int]]:
    len_seq = len(seq)
    if not len_seq:
        return "", []

    min_length_abs = round(len_seq*min_length)

    filter_list = [True] * len_seq
    i = 0

    while i < len_seq:
        if qscore[i] < threshold_q:
            filter_list[i] = False
            j = i + 1
            window_sum = qscore[i]
            window_size = 1

            while j < len_seq and (window_sum / window_size) < threshold_q:
                window_sum += qscore[j]
                window_size += 1
                filter_list[j] = False
                j += 1

            i = j
        else:
            i += 1

    if sum(filter_list) >= min_length_abs:
        return (''.join(s for s, f in zip(seq, filter_list) if f),
                [q for q, f in zip(qscore, filter_list) if f])
    return "", []


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

def fn_comp_rev(nameseq):
    nameseq = nameseq.replace('\n','')
    table= str.maketrans('ACTG','TGAC')
    comp = nameseq.translate(table)
    comp_rev= nameseq.translate(table)[::-1]
    return comp,comp_rev

def divide_into_kmer(sequence,k):
    kmer_list = []
    for i in range(0, len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmer_list.append(kmer)
    return kmer_list

def load_table(csv_path):
    df = pd.read_csv(csv_path)
    return df




