import argparse


import pandas as pd

from Sequence import Sequence
import ErroriPersonalizzati
from multiprocessing import Pool,cpu_count

"""
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', help='File Query di riferimento', required=True)
    parser.add_argument('-s', '--subject', help='File Subject di riferimento', required=True)
    parser.add_argument('-k','--kmer-length',type=int,default=22,help='Dimensione del kmer da estrarre')
    return parser.parse_args()

"""
query = Sequence('query.fasta')
list_partenza_query =query.parse_file()
kmer_query_list = query.kmer_indexing(22)
kmer_comprev_query_list = query.kmer_indexing_comp_rev(22)



query_1 = kmer_query_list[0:2]
query2 = kmer_query_list[2:]
query1r = kmer_comprev_query_list[0:2]
query2r = kmer_comprev_query_list[2:]
#print(query1,query2)


if kmer_query_list is None and kmer_comprev_query_list is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_query_list,list) or not isinstance(kmer_comprev_query_list,list):
    raise ErroriPersonalizzati.NotADict()


sub = Sequence('ref.fa')
list_partenza_subject = sub.parse_file()
kmer_subject_list = sub.kmer_indexing(22)
#print(kmer_subject_list)
kmer_comprev_subject_list = sub.kmer_indexing_comp_rev(22)


def create_query_df(kmer_query_list):
    headers = []
    kmers = []
    for element in range(0,len(kmer_query_list)):
        if element == 0 or element % 2 == 0:
            headers.append(kmer_query_list[element])
        else:
            kmers.append(kmer_query_list[element])
    df = pd.DataFrame('',index=kmers,columns=headers)
    df = df.reset_index()
    df.columns.values[0] = 'kmer'
    return df

def fill_df_query(df,nome_colonna):
    df[nome_colonna] = range(len(df))
    return df

def create_sub_df(kmer_sub_list):
    headers = []
    kmers = set()
    for element in range(0,len(kmer_sub_list)):
        if element == 0 or element % 2 == 0:
            headers.append(kmer_sub_list[element])
        else:
            kmer_sub_set = kmer_sub_list[element]
            kmers.update(kmer_sub_set)
    sorted_kmers = sorted(kmers)
    df = pd.DataFrame('',index=sorted_kmers,columns=headers)
    df = df.reset_index()
    df.columns.values[0] = 'kmer'
    return df
'''
def fill_sub_df(df,kmer_subject_list):
    for i in range(0,len(kmer_subject_list),2):
        coppia = kmer_subject_list[i:i+2]
        header = coppia[0]
        kmer_list = coppia[1]
        if header in df.columns[1:].tolist():
            for kmer in df['kmer']:
                if kmer in kmer_list:
                    if kmer_list.count(kmer) == 1:
                        df.at[kmer,header] = kmer_list.index(kmer)
                    elif kmer_list.count(kmer) > 1:
                        raise ValueError("aaaaa") # da rivedere
                else:
                    df.at[kmer,header] = None
    return df
    '''


def process_header(args):
    """
    Elaborazione per un singolo header.
    """
    header, kmer_list, df = args
    if header in df.columns[1:].tolist():  # Controlla se l'header è una colonna del DataFrame
        for kmer in df['kmer']:
            if kmer in kmer_list:
                if kmer_list.count(kmer) == 1:
                    df.at[kmer, header] = kmer_list.index(kmer)
                else:
                    raise ValueError("Il k-mer si ripete più di una volta in kmer_list.")
            else:
                df.at[kmer, header] = None
    return df[header]


# Funzione principale per parallelizzare il lavoro
def parallel_fill_sub_df(df, kmer_subject_list):
    """
    Parallelizza l'elaborazione del riempimento delle celle.
    """
    # Crea una lista di argomenti per ogni header
    tasks = []
    for i in range(0, len(kmer_subject_list), 2):
        coppia = kmer_subject_list[i:i + 2]
        header = coppia[0]
        kmer_list = coppia[1]
        tasks.append((header, kmer_list, df.copy()))

    # Usa il pool di processi per parallelizzare il lavoro
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(process_header, tasks)

    # Combina i risultati delle colonne elaborate
    for i, header in enumerate(df.columns[1:]):  # Escludi la colonna 'kmer'
        df[header] = results[i]

    return df



dataf1 = create_query_df(query_1)
dataf2 = create_query_df(query2)
datafr1 = create_query_df(query1r)
datafr2 = create_query_df(query2r)

dataf1_fill = fill_df_query(dataf1,"b6635d67cb594473ddba9f8cfba5d13d")
dataf2_fill = fill_df_query(dataf2,"4516aa60a483dd8c7bbc57098c45f1a5")
datafr1_fill = fill_df_query(datafr1,"b6635d67cb594473ddba9f8cfba5d13d")
datafr2_fill = fill_df_query(datafr2,"4516aa60a483dd8c7bbc57098c45f1a5")

datasub = create_sub_df(kmer_subject_list)
#c = fill_sub_df(datasub,kmer_subject_list)

#filled_sub_df = parallel_fill_sub_df(datasub,kmer_subject_list)

def main():
    query = Sequence('query.fasta')
    list_partenza_query = query.parse_file()
    kmer_query_list = query.kmer_indexing(22)
    kmer_comprev_query_list = query.kmer_indexing_comp_rev(22)

    query_1 = kmer_query_list[0:2]
    query2 = kmer_query_list[2:]
    query1r = kmer_comprev_query_list[0:2]
    query2r = kmer_comprev_query_list[2:]
    # print(query1,query2)

    if kmer_query_list is None and kmer_comprev_query_list is None:
        raise ErroriPersonalizzati.EmptyDict()
    if not isinstance(kmer_query_list, list) or not isinstance(kmer_comprev_query_list, list):
        raise ErroriPersonalizzati.NotADict()

    sub = Sequence('ref.fa')
    list_partenza_subject = sub.parse_file()
    kmer_subject_list = sub.kmer_indexing(22)
    # print(kmer_subject_list)
    kmer_comprev_subject_list = sub.kmer_indexing_comp_rev(22)

    dataf1 = create_query_df(query_1)
    dataf2 = create_query_df(query2)
    datafr1 = create_query_df(query1r)
    datafr2 = create_query_df(query2r)

    dataf1_fill = fill_df_query(dataf1, "b6635d67cb594473ddba9f8cfba5d13d")
    dataf2_fill = fill_df_query(dataf2, "4516aa60a483dd8c7bbc57098c45f1a5")
    datafr1_fill = fill_df_query(datafr1, "b6635d67cb594473ddba9f8cfba5d13d")
    datafr2_fill = fill_df_query(datafr2, "4516aa60a483dd8c7bbc57098c45f1a5")

    datasub = create_sub_df(kmer_subject_list)
    filled_sub_df = parallel_fill_sub_df(datasub, kmer_subject_list)
if __name__ == '__main__':
    main()






#final_df = filled_sub_df.merge(dataf1_fill,on='kmer',how='on')





'''
if kmer_subject_list is None and kmer_comprev_subject_list is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_subject_list,list) or not isinstance(kmer_comprev_subject_list,list):
    raise ErroriPersonalizzati.NotADict()


def find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)->dict:
    """
        Identifies seeds between query and subject dictionaries and organizes their positions.

        This function takes dictionaries representing k-mers from query sequences and subject sequences (including their
        complementary reverse sequences) and identifies shared k-mers. It organizes the positions of these shared k-mers
        in both the query and subject sequences into a nested dictionary structure.

        Parameters:
        kmer_query_dict : dict
            A dictionary where keys represent query sequence identifiers, and values are dictionaries of k-mers mapped to their positions.

        kmer_subject_dict : dict
            A dictionary where keys represent subject sequence identifiers, and values are dictionaries of k-mers mapped to their positions.

        kmer_comprev_subject_dict : dict
            A dictionary where keys represent subject sequence identifiers, and values are dictionaries of reverse-complement k-mers mapped to their positions.

        Returns:
        dict
            A dictionary containing seeds as keys. The values are nested dictionaries organizing the positions
            of these k-mers in query and subject sequences.
        """

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
    '''
'''
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
                                        
    return seed_dict
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

'''



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
'''


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

