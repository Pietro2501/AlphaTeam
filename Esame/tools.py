import ErroriPersonalizzati

def calculate_score(kmer1: str, kmer2: str,scoring_matrix) -> int:
    if len(kmer1) != len(kmer2):
        raise ValueError("I due kmer devono avere la stessa lunghezza!")

    score = 0
    for i in range(len(kmer1)):
        base1, base2 = kmer1[i], kmer2[i]
        score += scoring_matrix[(base1,base2)]
    return score

def generate_words(kmer,scoring_matrix,threshold = 0,max_words = None):
    nucleotidi = ['A', 'C', 'G', 'T']
    words = set()

    score = calculate_score(kmer,kmer,scoring_matrix)
    print(f" Original kmer={kmer}, base_score={score}")
    if score >= threshold:
        words.add(score)

    for i in range(len(kmer)):
        original_base = kmer[i]
        for base in nucleotidi:
            if base == original_base:
                continue  #non devo fare nessuna sostituzione
            new_kmer = kmer[:i] + base + kmer[i+1:]


            #devo ora calcolare il nuovo punteggio
            new_score = calculate_score(kmer,new_kmer,scoring_matrix)
            print(f"[DEBUG] variant={new_kmer}, variant_score={new_score}")
            if new_score >= threshold:
                words.add(new_kmer)

            if max_words is not None and len(words) >= max_words:
                return list(words)

    return list(words)

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
    for i in range(0, len(sequence)-k+1,k-2):
        kmer = sequence[i:i+k]
        kmer_diz.setdefault(kmer, [])
        kmer_diz[kmer].append(i)
    return kmer_diz
