import ErroriPersonalizzati
import pandas as pd
def calculate_score(kmer1: str, kmer2: str,scoring_matrix) -> int:
    if len(kmer1) != len(kmer2):
        raise ValueError("I due kmer devono avere la stessa lunghezza!")
    score = 0
    for i in range(len(kmer1)):
        base1, base2 = kmer1[i], kmer2[i]
        score += scoring_matrix[(base1,base2)]
    return score

def divide_into_kmer(sequence,k):
    if k <= 0:
        raise ErroriPersonalizzati.KmerError()
    if k > len(sequence):
        raise ErroriPersonalizzati.KmerTooLong()
    kmer_list = []
    for i in range(0, len(sequence)-k+1):
        kmer = sequence[i:i+k]
        #kmer_diz.setdefault(kmer, [])
        kmer_list.append(kmer)
    return kmer_list

def load_table(csv_path):
    df = pd.read_csv(csv_path)
    return df

