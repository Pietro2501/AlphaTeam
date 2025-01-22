import argparse
import time
import pandas as pd
import tools
from Sequence import Sequence
import ErroriPersonalizzati



query = Sequence('query.fasta')
list_partenza_query =query.parse_file()
kmer_query_list = query.kmer_indexing(22)
kmer_comprev_query_list = query.kmer_indexing_comp_rev(22)

query_partenza = list_partenza_query[0]





query_1 = kmer_query_list[0:2]
query2 = kmer_query_list[2:]
query1r = kmer_comprev_query_list[0:2]
query2r = kmer_comprev_query_list[2:]



if kmer_query_list is None and kmer_comprev_query_list is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_query_list,list) or not isinstance(kmer_comprev_query_list,list):
    raise ErroriPersonalizzati.NotADict()


sub = Sequence('ref.fa')
list_partenza_subject = sub.parse_file()
#print(list_partenza_subject)
kmer_subject_list = sub.kmer_indexing(22)
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

def fill_sub_df(df,kmer_subject_list):
    total_header = len(kmer_subject_list) // 2
    completed_headers = 0
    for i in range(0,len(kmer_subject_list),2):
        coppia = kmer_subject_list[i:i+2]
        header = coppia[0]
        kmer_list = coppia[1]
        print(f"Inizio elaborazione: {completed_headers + 1} su {total_header} - Header: {header}")
        start_time = time.time()
        if header in df.columns[1:].tolist():
            for kmer in df['kmer'].tolist():
                if kmer in kmer_list:
                    if kmer_list.count(kmer) == 1:
                        #print(kmer_list.index(kmer))
                        index = kmer_list.index(kmer)
                        df.loc[df['kmer'] == kmer,header] = index
                    elif kmer_list.count(kmer) > 1:
                        raise ValueError("aaaaa") # da rivedere
                else:
                    df.loc[df['kmer'] == kmer,header] = None
        completed_headers += 1
        tempo = time.time() - start_time
        print(f"Completamento: {completed_headers} su {total_header} - Tempo Impiegato {tempo:.2f} secondi")
    return df
def create_df_with_positions(query_df,subject_df,filename):
    result_df = pd.merge(query_df,subject_df,on='kmer',how='inner')
    result = {}
    query_position_col = [col for col in query_df.columns if col != "kmer"][0]
    for index, row in result_df.iterrows():
        kmer = row["kmer"]
        query_pos = int(row[query_position_col])
        for col in subject_df.columns[1:]:
            if not pd.isna(row[col]):
                if kmer not in result:
                    result[kmer] = {}
                result[kmer][col] = (query_pos, row[col])
    final_df = pd.DataFrame.from_dict(result,orient='index')
    if 'kmer' in final_df.columns:
        final_df = final_df.drop(columns=['kmer'])
    final_df.to_csv(filename +'.csv')
    return final_df

########## CREAZIONE QUERY DF #############
dataf1 = create_query_df(query_1)
dataf1.to_csv('a.csv')
dataf2 = create_query_df(query2)
datafr1 = create_query_df(query1r)
datafr2 = create_query_df(query2r)

############## FASE DI FILL DEI QUERY DF ################
dataf1_fill = fill_df_query(dataf1,"b6635d67cb594473ddba9f8cfba5d13d")
dataf1_fill.to_csv('data1fill.csv')
dataf2_fill = fill_df_query(dataf2,"4516aa60a483dd8c7bbc57098c45f1a5")
datafr1_fill = fill_df_query(datafr1,"b6635d67cb594473ddba9f8cfba5d13d")
datafr2_fill = fill_df_query(datafr2,"4516aa60a483dd8c7bbc57098c45f1a5")

########## CREAZIONE SUB DF #############
datasub = create_sub_df(kmer_subject_list)

############## CARICAMENTO DF SUB FILLATO
data = tools.load_table('filled_sub_df.csv')

################ TABELLE PER RICERCA SEED #########################
df_final_1 = create_df_with_positions(dataf1_fill,data,'final1_df')

df_final_2 = create_df_with_positions(dataf2_fill,data,'final2_df')

def find_seeds(df,filename,kmer_length=22):

    df = df.map(lambda x: x if isinstance(x, tuple) else None)
    seed_results = []

    for col in df.columns:
        column_seeds = []
        col_values = df[col].tolist()
        n = len(col_values)

        for i in range(n):
            current_value = col_values[i]
            prev_value = col_values[i - 1] if i > 0 else None
            next_value = col_values[i + 1] if i < n - 1 else None

            if current_value is not None:
                query_start, subject_start = current_value
                query_end = query_start + kmer_length
                subject_end = subject_start + kmer_length

                # Aggiungi una lista uniforme con start e end
                if prev_value is None and next_value is not None:
                    column_seeds.append([query_start, subject_start, None, None])
                elif prev_value is not None and next_value is None:
                    if column_seeds and column_seeds[-1][2:] == [None, None]:
                        column_seeds[-1][2:] = [query_end, subject_end]
                elif prev_value is None and next_value is None:
                    column_seeds.append([query_start, subject_start, query_end, subject_end])

        seed_results.append(column_seeds)


    max_rows = max(len(seeds) for seeds in seed_results)
    new_df_data = {
        df.columns[i]: seed_results[i] + [None] * (max_rows - len(seed_results[i]))
        for i in range(len(df.columns))
    }
    new_df = pd.DataFrame(new_df_data)

    new_df.index = range(1, len(new_df) + 1)
    new_df.to_csv(filename +'.csv')
    return new_df
seeds1 = find_seeds(df_final_1,'query1seeds')
seeds2 = find_seeds(df_final_2,'query2seeds')


def find_seed(kmer_query_dict, kmer_subject_dict, kmer_comprev_subject_dict) -> dict:
    seed_dict = {}
    for key1, inner_dict in kmer_query_dict.items():
        for kmer1, pos1 in inner_dict.items():
            for key2, sub_dict in kmer_subject_dict.items():
                for kmer2, pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_dict:
                            seed_dict[kmer1] = {'query': {key1: pos1}, 'subject': {key2: pos2}}
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
transizione = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
trasversione = {
    'A': ['C', 'T'],
    'C': ['A', 'G'],
    'G': ['C', 'T'],
    'T': ['A', 'G']
}

def get_sequence(header, list_partenza_subject):
    """
    Recupera la sequenza corrispondente a un dato header.

    Parameters:
    header (str): L'header del subject.
    list_partenza_subject (list of tuples): Lista di tuple (header_subject, sequenza_subject).

    Returns:
    list: Lista di caratteri della sequenza, oppure None se non trovata.
    """
    for item in list_partenza_subject:
        if item[0] == header:
            return list(item[1])  # Converti la stringa in lista di caratteri
    print(f"Errore: Subject '{header}' non trovato in list_partenza_subject.")
    return None
def handle_gaps(finestra_mismatch, max_gap):
    score_with_gap = -max_gap
    sequences = [finestra_mismatch]
    return score_with_gap, sequences
def extend_seed_right(sequence_query_ext, sequence_sub_ext, k, x_max):
    score = 0
    mismatch_consecutivi = 0
    a = 0  # Indice per l'estensione

    for a in range(len(sequence_query_ext)):
        if a >= len(sequence_sub_ext):
            break  # Evita di superare la lunghezza delle sequenze

        if sequence_query_ext[a] == sequence_sub_ext[a]:
            mismatch_consecutivi = 0
            score += 1
        else:
            mismatch_consecutivi += 1
            chiave = sequence_query_ext[a]
            if chiave in transizione and transizione[chiave] == sequence_sub_ext[a]:
                score -= 1
            elif chiave in trasversione and sequence_sub_ext[a] in trasversione[chiave]:
                score -= 1
            if mismatch_consecutivi == x_max:
                print(f"Interruzione per mismatch consecutivi a posizione {a}")
                break

    # Calcola la finestra di mismatch
    finestra_mismatch = sequence_query_ext[max(a - (x_max - 1), 0):a + 1]
    print(f"Finestra mismatch: {finestra_mismatch}")
    score_with_gap, sequences = handle_gaps(finestra_mismatch, max_gap=3)

    if score_with_gap > -x_max:
        new_extension = sequence_query_ext[:a - (x_max - 1)] + list(sequences[0])
        extension_right = new_extension
        print(f"Estensione con gap: {extension_right}")
    else:
        extension_right = sequence_query_ext[:a - (x_max - 1)] if a >= (x_max - 1) else sequence_query_ext[:a + 1]
        print(f"Estensione senza gap: {extension_right}")

    return extension_right, score


def extend_seed(df_seeds, query_partenza, list_partenza_subject, kmer_length=22, x_max=3,):

    query_partenza_header = query_partenza[0]
    query_partenza_sequence = query_partenza[1]
    query_header, sequence_query = query_partenza_header, query_partenza_sequence
    sequence_query = list(sequence_query)  # Converti la stringa in lista di caratteri

    extended_seeds_dict = {}

    # Itera attraverso ogni colonna del DataFrame
    for col in df_seeds.columns:
        seeds = df_seeds[col].dropna().tolist()  # Ottieni la lista di seed, escludendo None
        extended_seeds = []

        # Recupera la sequenza del subject corrente utilizzando la funzione get_sequence
        sequence_sub = get_sequence(col, list_partenza_subject)
        if sequence_sub is None:
            # Se la sequenza non è trovata, salta questa colonna
            extended_seeds_dict[col] = extended_seeds
            continue

        # Itera attraverso ogni seed nella colonna
        for i in range(len(seeds)):
            current_seed = seeds[i]
            if not current_seed:
                continue  # Salta seed vuoti

            start_q, start_s, end_q, end_s = current_seed

            # Assicurati che start_q, start_s, end_q, end_s siano interi
            try:
                start_q = int(start_q)
                start_s = int(start_s)
                end_q = int(end_q)
                end_s = int(end_s)
            except ValueError:
                print(f"Errore: I valori di start_q, start_s, end_q, end_s devono essere interi. Seed: {current_seed}")
                continue

            # Controlla se c'è un seed successivo
            if i < len(seeds) - 1:
                next_seed = seeds[i + 1]
                if next_seed:
                    next_start_q, next_start_s, next_end_q, next_end_s = next_seed

                    # Assicurati che anche il next_seed abbia valori interi
                    try:
                        next_start_q = int(next_start_q)
                        next_start_s = int(next_start_s)
                        next_end_q = int(next_end_q)
                        next_end_s = int(next_end_s)
                    except ValueError:
                        print(f"Errore: I valori di next_start_q, next_start_s, next_end_q, next_end_s devono essere interi. Seed: {next_seed}")
                        continue

                    # Calcola le differenze
                    diff_q = next_start_q - end_q
                    diff_s = next_start_s - end_s

                    if diff_q == diff_s:
                        # Nessun gap, estendi il seed
                        sequence_query_ext = sequence_query[end_q:]
                        sequence_sub_ext = sequence_sub[end_s:]
                        extension_right, score = extend_seed_right(
                            sequence_query_ext,
                            sequence_sub_ext,
                            kmer_length,
                            x_max
                        )

                        # Aggiorna end_q e end_s
                        new_end_q = end_q + len(extension_right)
                        new_end_s = end_s + len(extension_right)

                        extended_seed = [start_q, start_s, new_end_q, new_end_s]
                        extended_seeds.append(extended_seed)
                    else:
                        # Gestisci i gap
                        if diff_q < diff_s:
                            # Gap nella query
                            gap_size = diff_s - diff_q
                            new_end_q = end_q + gap_size
                            new_end_s = end_s + gap_size

                            # Aggiungi il seed esteso con gap nella query
                            extended_seed = [start_q, start_s, new_end_q, end_s]
                            extended_seeds.append(extended_seed)

                            # Valuta l'estensione con una finestra
                            finestra_mismatch = sequence_query[new_end_q:new_end_q + gap_size]
                            score_with_gap, sequences = handle_gaps(finestra_mismatch, max_gap=3)
                            print(f"Gap nella query: {finestra_mismatch}, Score with gap: {score_with_gap}, Sequences: {sequences}")
                            # Puoi usare 'sequences' per ulteriori elaborazioni

                        else:
                            # Gap nel subject
                            gap_size = diff_q - diff_s
                            new_end_q = end_q + gap_size
                            new_end_s = end_s + gap_size

                            # Aggiungi il seed esteso con gap nel subject
                            extended_seed = [start_q, start_s, end_q, new_end_s]
                            extended_seeds.append(extended_seed)

                            # Valuta l'estensione con una finestra
                            finestra_mismatch = sequence_sub[new_end_s:new_end_s + gap_size]
                            score_with_gap, sequences = handle_gaps(finestra_mismatch, max_gap=3)
                            print(f"Gap nel subject: {finestra_mismatch}, Score with gap: {score_with_gap}, Sequences: {sequences}")
                            # Puoi usare 'sequences' per ulteriori elaborazioni
            else:
                # Seed finale, potrebbe non avere un successivo
                sequence_query_ext = sequence_query[end_q + kmer_length:]
                sequence_sub_ext = sequence_sub[end_s + kmer_length:]
                extension_right, score = extend_seed_right(
                    sequence_query_ext,
                    sequence_sub_ext,
                    kmer_length,
                    x_max
                )

                new_end_q = end_q + len(extension_right)
                new_end_s = end_s + len(extension_right)

                extended_seed = [start_q, start_s, new_end_q, new_end_s]
                extended_seeds.append(extended_seed)

        # Aggiungi i seed estesi alla lista del dizionario
        extended_seeds_dict[col] = extended_seeds

    # Converti il dizionario in un nuovo DataFrame
    # Trova la lunghezza massima tra le colonne
    max_rows = max(len(seeds) for seeds in extended_seeds_dict.values())

    # Uniforma la lunghezza di tutte le colonne aggiungendo righe vuote con un singolo None
    new_df_data = {
        col: extended_seeds_dict[col] + [None] * (max_rows - len(extended_seeds_dict[col]))
        for col in df_seeds.columns
    }

    new_df = pd.DataFrame(new_df_data)

    # Imposta l'indice del nuovo DataFrame
    new_df.index = range(1, len(new_df) + 1)

    return new_df

extended_seeds = extend_seed(seeds1,query_partenza,list_partenza_subject,22,6)







#p = extend_seed(seeds1,list_partenza_query,list_partenza_subject)
#final_1
#(379,52)
# final_2
# (381,49)



