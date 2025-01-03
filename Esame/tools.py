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
