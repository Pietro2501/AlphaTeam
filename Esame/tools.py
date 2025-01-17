import ErroriPersonalizzati

def calculate_score(kmer1: str, kmer2: str,scoring_matrix) -> int:
    if len(kmer1) != len(kmer2):
        raise ValueError("I due kmer devono avere la stessa lunghezza!")
    score = 0
    for i in range(len(kmer1)):
        base1, base2 = kmer1[i], kmer2[i]
        score += scoring_matrix[(base1,base2)]
    return score



def divide_into_kmer(sequence,k):
    """
    Divide la sequenza in kmer di lunghezza k, restituendo un dizionario
    dove la chiave sono i kmer e il valore è una lista di posizioni.

    Verranno sollevati errori personalizzati se:
    - k è minore o uguale a zero
    - k è maggiore della lunghezza della sequenza

    Parametri:
    sequence: str
        Sequenza di caratteri
    k: int
        Dimensione fissa del kmer da estrarre a partire dalla sequenza

    Return:
    kmer_diz: dict
        Un dizionario in cui le chiavi sono i kmer e i valori sono
        le liste di posizioni.
    """
    if k <= 0:
        raise ErroriPersonalizzati.KmerError()
    if k > len(sequence):
        raise ErroriPersonalizzati.KmerTooLong()
    kmer_diz = {}
    for i in range(0, len(sequence)-k+1,k):
        kmer = sequence[i:i+k]
        kmer_diz.setdefault(kmer, [])
        kmer_diz[kmer].append(i)
    return kmer_diz
