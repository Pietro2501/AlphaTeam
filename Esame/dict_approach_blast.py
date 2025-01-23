'''
if kmer_subject_list is None and kmer_comprev_subject_list is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_subject_list, list) or not isinstance(kmer_comprev_subject_list, list):
    raise ErroriPersonalizzati.NotADict()
'''
'''

"""
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', help='File Query di riferimento', required=True)
    parser.add_argument('-s', '--subject', help='File Subject di riferimento', required=True)
    parser.add_argument('-k','--kmer-length',type=int,default=22,help='Dimensione del kmer da estrarre')
    return parser.parse_args()

"""
def find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)->dict:
    

    seed_dict = {}
    for key1,inner_dict in kmer_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_dict:
                            seed_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_dict[kmer1]['query']:
                                seed_dict[kmer1]['query'][key1] = pos1
                            else:
                                for p in pos1:
                                        if p not in seed_dict[kmer1]['query'][key1]:
                                            seed_dict[kmer1]['query'][key1].append(p)

                            if key2 not in seed_dict[kmer1]['subject']:
                                seed_dict[kmer1]['subject'][key2] = pos2
                            else:
                                for p in pos2:
                                    if p not in seed_dict[kmer1]['subject'][key2]:
                                        seed_dict[kmer1]['subject'][key2].append(p)

    for key1,inner_dict in kmer_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_comprev_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_dict:
                            seed_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_dict[kmer1]['query']:
                                seed_dict[kmer1]['query'][key1] = pos1
                            else:
                                for p in pos1:
                                        if p not in seed_dict[kmer1]['query'][key1]:
                                            seed_dict[kmer1]['query'][key1].append(p)

                            if key2 not in seed_dict[kmer1]['subject']:
                                seed_dict[kmer1]['subject'][key2] = pos2
                            else:
                                for p in pos2:
                                    if p not in seed_dict[kmer1]['subject'][key2]:
                                        seed_dict[kmer1]['subject'][key2].append(p)
                                        '''
'''
def build_seed_list_from_seed_dict(seed_dict):
    """
    Trasforma il dizionario dei seed (che ha come chiavi i kmer, e come valori
    i dizionari di posizioni in query e subject) in una lista di tuple
    (query_id, subject_id, start_query, start_subject, k).

    Esempio di seed_dict:
      {
        'TGAGGAATATTGGTCAATGGGC': {
           'query': {
               'b6635d67cb594473ddba9f8cfba5d13d': [0, 5, ...]
            },
           'subject': {
               'ref_id_1': [100, 310, ...],
               'ref_id_2': [45, ...]
           }
        },
        ...
      }

    Output:
      [
        ('b6635d67cb594473ddba9f8cfba5d13d', 'ref_id_1', 0, 100, lunghezza_kmer),
        ('b6635d67cb594473ddba9f8cfba5d13d', 'ref_id_1', 5, 310, lunghezza_kmer),
        ...
      ]
    """
    seed_list = []

    # Scorri tutti i kmer
    for kmer, info_dict in seed_dict.items():
        # lunghezza del kmer
        k = len(kmer)

        # info_dict dovrebbe essere fatto così:
        # {
        #   'query':   { qID: [posQ1, posQ2, ...], ... },
        #   'subject': { sID: [posS1, posS2, ...], ... }
        # }
        query_part = info_dict.get('query', {})
        subject_part = info_dict.get('subject', {})

        # Scorri tutti i queryID e le loro posizioni
        for qID, qPositions in query_part.items():
            for posQ in qPositions:
                # Scorri tutti i subjectID e le loro posizioni
                for sID, sPositions in subject_part.items():
                    for posS in sPositions:
                        # Crea la tupla finale per ciascuna combinazione
                        seed_tuple = (qID, sID, posQ, posS, k)
                        seed_list.append(seed_tuple)

    return seed_list
'''
'''
a = find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)
#print(a)

seed_list = build_seed_list_from_seed_dict(a)
print(seed_list)

schema = []
for kmer, inner_dict in a.items():
    for que_sub, sub_dict in inner_dict.items():
        if que_sub == 'query':
            for j in sub_dict.items():
                query = j
        if que_sub == 'subject':
            for i in sub_dict.items():
                subject = i
                schema.append(kmer)
                schema.append(query)
                schema.append(subject)

def find_consecutive_seeds(schema):
    for i in range(0,len(schema),3):
        prova = schema[i:i + 3]
        query = prova[1][0]
        subject = prova[2][0]
        pos_query = prova[1][1]
        pos_sub = prova[2][1]
        kmer = prova[0]




""" #CHIEDERE AL PROF SE è NECESSARIO
def find_comprev_seed(kmer_comprev_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)->dict:
    seed_comprev_dict = {}
    for key1,inner_dict in kmer_comprev_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_comprev_dict:
                            seed_comprev_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_comprev_dict[kmer1]['query']:
                                seed_comprev_dict[kmer1]['query'][key1] = [pos1]
                            if key2 not in seed_comprev_dict[kmer1]['subject']:
                                seed_comprev_dict[kmer1]['subject'][key2] = [pos2]

    for key1,inner_dict in kmer_comprev_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_comprev_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_comprev_dict:
                            seed_comprev_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_comprev_dict[kmer1]['query']:
                                seed_comprev_dict[kmer1]['query'][key1] = [pos1]
                            if key2 not in seed_comprev_dict[kmer1]['subject']:
                                seed_comprev_dict[kmer1]['subject'][key2] = [pos2]

    return seed_comprev_dict

"""


s = 38
x = 6
schema = []
for kmer, inner_dict in a.items():
    for que_sub, sub_dict in inner_dict.items():
        if que_sub == 'query':
            for j in sub_dict.items():
                query = j
        if que_sub == 'subject':
            for i in sub_dict.items():
                subject = i
                schema.append(kmer)
                schema.append(query)
                schema.append(subject)

#print(schema)


transizione = {'A':'G','G':'A','C':'T','T':'C'}
trasversione = {
    'A': ['C', 'T'],
    'C': ['A', 'G'],
    'G': ['C', 'T'],
    'T': ['A', 'G']
}

def extend_seed_right(sequence_query, sequence_sub, start_query, start_subject,k,x_max):
    """
        Extends a seed alignment to the right between query and subject sequences.

        This function extends a seed match between a query sequence and a subject sequence
        to the right, while scoring the alignment. It stops the extension if a maximum number
        of consecutive mismatches (`x_max`) is reached.

        Parameters:
        sequence_query : str
            The query sequence in which the seed alignment starts.

        sequence_sub : str
            The subject sequence in which the seed alignment starts.

        start_query : int
            The starting position of the seed in the query sequence.

        start_subject : int
            The starting position of the seed in the subject sequence.

        k : int
            The length of the initial seed.

        x_max : int
            The maximum number of consecutive mismatches allowed before terminating the extension.

        Returns:
        tuple
            A tuple containing:
            - extension_right (str): The extended matching sequence to the right.
            - score (int): The alignment score based on matches, transitions, and transversions.
        """
    score = 0
    mismatch_consecutivi = 0

    sequence_query = sequence_query[start_query + k:]
    sequence_sub = sequence_sub[start_subject + k:]

    for a in range(0, len(sequence_query)):
        if sequence_query[a] == sequence_sub[a]:
            mismatch_consecutivi = 0
            score += 1
        else:
            mismatch_consecutivi += 1
            chiave = sequence_query[a]
            if transizione[chiave] == sequence_sub[a]:
                score -= 1
            elif trasversione[chiave] == sequence_sub[a]:
                score -= 1
            if mismatch_consecutivi == x_max:
                break

    finestra_mismatch = sequence_query[a - (x_max - 1):a + 1]
    pos_iniziale_finestra = finestra_mismatch[0]
    score_finestra = -x_max
    score_with_gap,sequences = handle_gaps(finestra_mismatch,3)

    if score_with_gap > score_finestra:
        new_extension = sequence_query.replace(finestra_mismatch,sequences[0])
        new_extension = new_extension[pos_iniziale_finestra:]


    extension_right = sequence_query[0:a-(x_max-1)]

    return extension_right,score

def calculate_score(sequence_query,sequence_sub):
    score = 0
    for b in range(len(sequence_sub)):
        if b >= len(sequence_query):
            score -= 2
            continue
        query_char = sequence_query[b]
        sub_char = sequence_sub[b]

        if query_char == sub_char:
            score += 1
        else:
            if query_char == '_':
                if b == 0 or sequence_query[b - 1] != '_':
                    score -= 2
                else:
                    score -= 1
            elif query_char in transizione and transizione[query_char] == sub_char:
                score -= 1
            elif query_char in trasversione and sub_char in trasversione[query_char]:
                score -= 1
            else:
                score -= 1
    return score


def handle_gaps(finestra_mismatch, soglia_gap):
    sequence_sub_finestra = finestra_mismatch
    gap_scores = []
    sequence_list = []

    for a in range(soglia_gap):
        sequence_query_finestra = '_' * (a + 1) + finestra_mismatch
        sequence_list.append(sequence_query_finestra)
        gap_score = calculate_score(sequence_query_finestra,sequence_sub_finestra)
        gap_scores.append(gap_score)


    max_score = max(gap_scores)
    max_indices = [i for i, score in enumerate(gap_scores) if score == max_score]
    sequences =[sequence_list[i] for i in max_indices]
    return max_score, sequences

print(handle_gaps('ACGTCG',3))

def extend_seed_left(sequence_query, sequence_sub, start_query, start_subject, k, x_max):
    """
        Extends a seed alignment to the left between query and subject sequences.

        This function extends a seed match to the left, starting from specified positions in the query
        and subject sequences. The extension stops when the maximum number of consecutive mismatches (`x_max`)
        is reached or when the beginning of either sequence is encountered. The alignment score is computed
        based on matches, transitions, and transversions.

        Parameters:
        sequence_query : str
            The query sequence in which the seed alignment starts.

        sequence_sub : str
            The subject sequence in which the seed alignment starts.

        start_query : int
            The starting position of the seed in the query sequence.

        start_subject : int
            The starting position of the seed in the subject sequence.

        k : int
            The length of the initial seed.

        x_max : int
            The maximum number of consecutive mismatches allowed before terminating the extension.

        Returns:
        tuple
            A tuple containing:
            - estensione_sinistra (str): The extended matching sequence to the left.
            - score (int): The alignment score based on matches, transitions, and transversions.

        """
    score = 0
    match_consecutivi = 0

    sequence_query_left = sequence_query[:start_query]
    sequence_subject_left = sequence_sub[:start_subject]


    estensione_sinistra = ""
    mismatch_consecutivi = 0

    if start_query == 0:
        estensione_sinistra = ""
        #print("Nessuna estensione possibile a sinistra")
    else:
        for queryBase, subjectBase in zip(
                reversed(sequence_query_left),
                reversed(sequence_subject_left)
        ):
            if queryBase == subjectBase:
                score += 1
                mismatch_consecutivi = 0
            else:
                mismatch_consecutivi += 1
                chiave = queryBase
                if transizione[chiave] == subjectBase:
                    score -= 1
                elif trasversione[chiave] == subjectBase:
                    score -= 1

                if mismatch_consecutivi == x_max:
                    #print(f"Estensione a sinistra:  Mi sono fermato dopo {mismatch_consecutivi} mismatch consecutivi.")
                    break


            estensione_sinistra = queryBase + estensione_sinistra

    return estensione_sinistra, score


def extend_seed(schema, diz_partenza_query, diz_partenza_subject):
    """
        Extends seed alignments based on a given schema and computes alignment scores and positions.

        This function takes a schema that defines seed alignments between query and subject sequences,
        extends the alignments both to the right and left, and calculates the resulting alignment scores.
        It returns a dictionary containing extended alignments, their positions, and scores.

        Parameters:
        schema : list
            A list of tuples representing the schema of seed alignments. Each tuple consists of:
            - hsp (str): The initial seed sequence.
            - query_info (tuple): A tuple with the query identifier and its seed positions ([start, end]).
            - subject_info (tuple): A tuple with the subject identifier and its seed positions ([start, end]).

        diz_partenza_query : dict
            A dictionary where keys are query identifiers and values are the corresponding query sequences.

        diz_partenza_subject : dict
            A dictionary where keys are subject identifiers and values are the corresponding subject sequences.

        Returns:
        dict
            A dictionary where keys are the extended high-scoring pairs (HSPs), and values are lists containing:
            - Query information: A dictionary with the query identifier and the start and end positions of the extended alignment.
            - Subject information: A dictionary with the subject identifier and the start and end positions of the extended alignment.
            - Score information: A dictionary with the total alignment score.

        """
    hsp_dict = {}
    for i in range(0,len(schema),3):

        prova = schema[i:i+3]
        query = prova[1][0]
        subject = prova[2][0]
        pos_query = prova[1][1]
        pos_sub = prova[2][1]
        hsp = prova[0]
        k = len(hsp)
        start_query = pos_query[0]
        start_subject = pos_sub[0]

        for j in diz_partenza_query.keys():
            if query == j:
                sequence_query = diz_partenza_query[j]

        for z in diz_partenza_subject.keys():
            if subject == z:
                sequence_sub = diz_partenza_subject[z]

        right_ext = extend_seed_right(sequence_query, sequence_sub, start_query, start_subject, k=22, x_max=6)[0]
        left_ext = extend_seed_left(sequence_query,sequence_sub,start_query,start_subject,k=22,x_max=6)[0]
        right_score = extend_seed_right(sequence_query, sequence_sub, start_query, start_subject, k=22, x_max=6)[1]
        left_score = extend_seed_left(sequence_query, sequence_sub, start_query, start_subject, k=22, x_max=6)[1]

        hsp = left_ext+hsp+right_ext
        score = left_score + k + right_score

        pos_init_hsp_query = start_query- len(left_ext)
        pos_init_hsp_sub = start_subject - len(left_ext)
        pos_fin_hsp_query = start_query + k + len(right_ext)
        pos_fin_hsp_sub = start_subject + k + len(right_ext)


        hsp_dict[hsp] = [{query:[pos_init_hsp_query,pos_fin_hsp_query]},{subject:[pos_init_hsp_sub,pos_fin_hsp_sub]},{'score totale': score}]

    return hsp_dict

#print(extend_seed(schema, diz_partenza_query, diz_partenza_subject))





'''
'''
schema = []
kmer = "TGAGGAATATTGGTCAATGGGC"
query_header = "b6635d67cb594473ddba9f8cfba5d13d"
subject_header = "MJ030-2-barcode67-umi101484bins-ubs-3"
pos_query = [0]
pos_subject = [330]

schema.append(kmer)
schema.append((query_header, pos_query))
schema.append((subject_header, pos_subject))

k = len(kmer)
x_max = 6


hsp_extended,score_extended = extend_seed_left(
    schema,
    diz_partenza_query_def,
    diz_partenza_subject_def,
    k,
    x_max
)
'''
'''
hsp_right,score_right = extend_seed_right(
    diz_partenza_query,
    diz_partenza_subject,
    pos_query,
    pos_subject,
    k,
    x_max
)
print("\nRisultato a sinistra")
print(kmer)
print(f"HSP Esteso senza seed: {hsp_extended}")
print(f"Score: {score_extended}")

print("\nRisultato a destra")
print(kmer)
print(f"HSP Esteso senza seed: {hsp_right}")
print(f"Score: {score_right}")

print("\nRisultato finale")
print(kmer)
print(hsp_extended+f'\033[1m{kmer}\033[0m'+hsp_right)
print(score_extended+score_right)










"""
def main():
    args = parse_args()

    query = Query(args.query)
    query.parse_file()

    subject = Subject(args.subject)
    subject.parse_file()

    kmer_query_dict = query.kmer_indexing(args.kmer_length)
    kmer_comprev_query_dict = query.kmer_indexing_comp_rev(args.kmer_length)

    if kmer_query_dict is None and kmer_comprev_query_dict is None:
        raise ErroriPersonalizzati.EmptyDict()
    if not isinstance(kmer_query_dict, dict) or not isinstance(kmer_comprev_query_dict, dict):
        raise ErroriPersonalizzati.NotADict()

    kmer_subject_dict = subject.subject_indexing(args.kmer_length)
    kmer_comprev_subject_dict = subject.subject_indexing_comp_rev(args.kmer_length)

    if kmer_subject_dict is None and kmer_comprev_subject_dict is None:
        raise ErroriPersonalizzati.EmptyDict()
    if not isinstance(kmer_subject_dict, dict) or not isinstance(kmer_comprev_subject_dict, dict):
        raise ErroriPersonalizzati.NotADict()

    a = find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)
    b = find_comprev_seed(kmer_comprev_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)

    print(f"Seed del forward: \033[1m{a}\033[0m")
    print()
    print(f"Seed del revertito complementare: \033[1m{b}\033[0m")

if __name__ == "__main__":
    main()
'''

import pandas as pd

def extend_seed(df_seeds, query_partenza, list_partenza_subject, kmer_length=22, x_max=3):
    
    query_header = query_partenza[0]
    sequence_query = query_partenza[1]
    
    extended_seeds_dict = {}
    scores_dict = {}

    # Itera attraverso ogni colonna del DataFrame
    for col in df_seeds.columns:
        seeds = df_seeds[col].dropna().tolist()  # Ottieni la lista di seed, escludendo None
        extended_seeds = []
        scores = []

        # Recupera la sequenza del subject corrente
        sequence_sub = get_sequence(col, list_partenza_subject)
        if sequence_sub is None:
            # Se la sequenza non è trovata, salta questa colonna
            extended_seeds_dict[col] = extended_seeds
            scores_dict[col] = scores
            continue

        if len(seeds) > 0:
            # Itera attraverso ogni seed nella colonna
            for i in range(len(seeds)):
                current_seed = seeds[i]
                if isinstance(current_seed, tuple):
                    start_q, start_s, end_q, end_s = map(int, current_seed)
                else:
                    start_q, start_s, end_q, end_s = map(int, current_seed.split(","))

                # Estensione iniziale
                hsp_query = sequence_query[start_q:end_q]
                hsp_sub = sequence_sub[start_s:end_s]

                # Se non ci sono altri seed successivi o siamo sull'ultimo seed
                if i == len(seeds) - 1:
                    extended_seeds.append((start_q, start_s, end_q, end_s))
                    scores.append(len(hsp_query))  # Usa la lunghezza come score (adatta secondo le tue esigenze)
                    continue

                # Recupera il prossimo seed
                next_seed = seeds[i + 1]
                if isinstance(next_seed, tuple):
                    next_start_q, next_start_s, next_end_q, next_end_s = map(int, next_seed)
                else:
                    next_start_q, next_start_s, next_end_q, next_end_s = map(int, next_seed.split(","))

                # Calcola le differenze tra la fine del seed corrente e l'inizio del successivo
                diff_q = next_start_q - end_q
                diff_s = next_start_s - end_s

                if diff_q == diff_s:
                    # Nessun gap, estendi direttamente
                    sequence_query_ext = sequence_query[end_q:next_start_q]
                    sequence_sub_ext = sequence_sub[end_s:next_start_s]

                else:
                    # Gap presente, determina dove inserirlo
                    gap_size = abs(diff_q - diff_s)

                    if diff_q < diff_s:
                        # Gap nella query
                        sequence_query_ext = sequence_query[end_q:]
                        sequence_sub_ext = sequence_sub[end_s:end_s + gap_size] + sequence_sub[end_s + gap_size:next_start_s]

                    else:
                        # Gap nel subject
                        sequence_query_ext = sequence_query[end_q:end_q + gap_size] + sequence_query[end_q + gap_size:next_start_q]
                        sequence_sub_ext = sequence_sub[end_s:]

                # Prova l'estensione verso destra
                extension_query = sequence_query[end_q:next_start_q]
                extension_sub = sequence_sub[end_s:next_start_s]

                if len(extension_query) == len(extension_sub):
                    hsp_query += extension_query
                    hsp_sub += extension_sub

                # Calcola il nuovo HSP
                new_end_q = end_q + len(extension_query)
                new_end_s = end_s + len(extension_sub)

                extended_seeds.append((start_q, start_s, new_end_q, new_end_s))
                scores.append(len(hsp_query))  # Usa la lunghezza come score

        extended_seeds_dict[col] = extended_seeds
        scores_dict[col] = scores

    # Creazione del DataFrame finale
    max_rows = max(len(extended_seeds_dict[col]) for col in extended_seeds_dict)

    # Unisci HSP e scores in un unico DataFrame
    new_df_data = {}
    for col in df_seeds.columns:
        extended_seeds = extended_seeds_dict.get(col, [])
        scores = scores_dict.get(col, [])

        # Allinea la lunghezza delle liste
        extended_seeds += [None] * (max_rows - len(extended_seeds))
        scores += [None] * (max_rows - len(scores))

        new_df_data[col] = extended_seeds
        new_df_data[f"{col}_score"] = scores

    new_df = pd.DataFrame(new_df_data)
    return new_df

# Funzioni placeholder per get_sequence e calculate_score (da implementare secondo il tuo contesto)
def get_sequence(header, list_subjects):
    for subject in list_subjects:
        if subject[0] == header:
            return subject[1]
    return None

def calculate_score(query, subject):
    # Calcola lo score come lunghezza dell'HSP (puoi sostituire con il tuo metodo di scoring)
    return len(query)

import pandas as pd

def test_extend_seed_realistic():
    print("Running realistic test cases for `extend_seed`...\n")
    
    # Caso 1: Seed semplice con estensione su regioni non identiche
    print("Caso 1: Seed semplice con estensione su regioni non identiche")
    df_case1 = pd.DataFrame({
        "subject_1": [["10,20,20,30"]],  # Seed iniziale
    })
    query_case1 = ["query_1", "ACGTTGACCTAGCGTAGCTAGCGTTAGCGTACGTA"]
    list_subjects_case1 = [["subject_1", "TGCATGACCTAGGGTAGCCCGGTAGCGTCCGTACG"]]
    result1 = extend_seed(df_case1, query_case1, list_subjects_case1)
    print(result1, "\n")
    
    # Caso 2: Due seed distanti con possibili gap su query e subject
    print("Caso 2: Due seed distanti con gap su query e subject")
    df_case2 = pd.DataFrame({
        "subject_1": [["5,15,10,20", "30,40,35,45"]],  # Seed con gap
    })
    query_case2 = ["query_2", "AGCTTAGCGTTAGCCCGTACCGTAGCTAGCTAGCTGCTG"]
    list_subjects_case2 = [["subject_1", "AGGTTAGCGTGGCCTCGTACCATTAGGTAGCTAGC"]]
    result2 = extend_seed(df_case2, query_case2, list_subjects_case2)
    print(result2, "\n")
    
    # Caso 3: Seed multipli, con estensioni sovrapposte
    print("Caso 3: Seed multipli, con estensioni sovrapposte")
    df_case3 = pd.DataFrame({
        "subject_1": [["10,20,20,30", "25,35,35,45"]],  # Overlapping seeds
    })
    query_case3 = ["query_3", "TGCATAGCGTTAGCTAGCCCTTACGATGCCGATAC"]
    list_subjects_case3 = [["subject_1", "TACGTAGCGTTCGCTAGCCCTTGGTGCCGATAGCC"]]
    result3 = extend_seed(df_case3, query_case3, list_subjects_case3)
    print(result3, "\n")
    
    # Caso 4: Seed con estensione parziale su regioni divergenti
    print("Caso 4: Seed con estensione parziale su regioni divergenti")
    df_case4 = pd.DataFrame({
        "subject_1": [["0,0,5,5"]],  # Estensione in regioni non perfette
    })
    query_case4 = ["query_4", "AACGTAGTTACGGGTACGATTAGCCTGACCTGAC"]
    list_subjects_case4 = [["subject_1", "AACGTACTTACGCGGACGCTTAGGATGACATGAC"]]
    result4 = extend_seed(df_case4, query_case4, list_subjects_case4)
    print(result4, "\n")
    
    # Caso 5: Nessun seed presente (controllo su input vuoti)
    print("Caso 5: Nessun seed presente (input vuoti)")
    df_case5 = pd.DataFrame({
        "subject_1": [[]],  # Nessun seed
    })
    query_case5 = ["query_5", "AGCTTAGCGTTAGCCCGTACCGTAGCTAGCTAGCTGCTG"]
    list_subjects_case5 = [["subject_1", "AGGTTAGCGTGGCCTCGTACCATTAGGTAGCTAGC"]]
    result5 = extend_seed(df_case5, query_case5, list_subjects_case5)
    print(result5)
  