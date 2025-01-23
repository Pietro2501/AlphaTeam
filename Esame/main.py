import argparse
import time
import pandas as pd
import tools
from Sequence import Sequence
import ErroriPersonalizzati
import os


query = Sequence('C:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam\\Esame\\query.fasta')
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


sub = Sequence('C:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam\\Esame\\ref.fa')
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
data = tools.load_table('C:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam\\Esame\\filled_sub_df.csv')

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
            return item[1] # Converti la stringa in lista di caratteri
    print(f"Errore: Subject '{header}' non trovato in list_partenza_subject.")
    return None

def calculate_score(sequence_1,sequence_2):
    score = 0
    stop_calc = min(len(sequence_1),len(sequence_2))
    flag = 0
    for b in range(stop_calc):
        query_char = sequence_1[b]
        sub_char = sequence_2[b]
        #print(b,query_char,sub_char)
        
        if query_char == sub_char:
            score += 1
            flag = 0
        else:
            if query_char == '_' or sub_char == '_':
                if flag == 0:
                    score -= 2
                    flag = 1
                else:
                    score -= 1
            elif query_char in transizione and transizione[query_char] == sub_char:
                score -= 1
                flag = 0
            elif query_char in trasversione and sub_char in trasversione[query_char]:
                score -= 1
                flag = 0
        #print(f'score= {score}')
    return score
  
def handle_gaps(finestra_mismatch_query, finestra_mismatch_sub, gap_size, gap_target, x_max):
    """
    La funzione gestisce i gap. Se sappiamo che sono sulla query, li aggiunge lì; se sono sulla subject, il contrario; altrimenti prova su entrambi.
    """
    max_gap = 3
    sequence_query_list = []
    sequence_sub_list = []
    gap_numbers = min(gap_size,max_gap)

    if gap_target == 'query':
        for a in range(gap_numbers):
            sequence_finestra_gap = '_' * (a + 1) + finestra_mismatch_query  # Gap a sinistra (query)
            gap_score = calculate_score(sequence_finestra_gap, finestra_mismatch_sub)
            if gap_score >= -x_max:
                sequence_query_list.append(sequence_finestra_gap)
            sequence_sub_list.append(finestra_mismatch_sub)

    if gap_target == 'subject':
        for a in range(gap_numbers):
            sequence_finestra_gap = '_' * (a + 1) + finestra_mismatch_sub   # Gap a destra (subject)
            gap_score = calculate_score(sequence_finestra_gap, finestra_mismatch_query)
            if gap_score >= -x_max:
                sequence_sub_list.append(sequence_finestra_gap)  
            sequence_query_list.append(finestra_mismatch_query)     

    return sequence_query_list, sequence_sub_list

def extend_seed_right(sequence_query_ext, sequence_sub_ext, gap_size,perform_gap,gap_target, x_max): ###FUNZIONA MA CERCHIAMO DI IMPLEMENTARLA
    mismatch_consecutivi = 0
    last_valid_index = 0
    #print(len(sequence_query_ext),len(sequence_sub_ext))
    for a in range(len(sequence_query_ext)):# a indice per l'estensione
        if a >= len(sequence_sub_ext):
            break  # Evita di superare la lunghezza delle sequenze
        #print(a,sequence_query_ext[a],sequence_sub_ext[a])
        if sequence_query_ext[a] == sequence_sub_ext[a]:
            mismatch_consecutivi = 0
        else:
            mismatch_consecutivi += 1
            
            if mismatch_consecutivi == x_max:
                #print(f"Interruzione per mismatch consecutivi a posizione {a}")
                last_valid_index = a
                break

    # Calcola la finestra di mismatch se query e subject hanno lunghezze differenti
    extension_right_query=[]
    extension_right_sub = []

    if perform_gap:
        if mismatch_consecutivi == x_max:
            finestra_mismatch_query = sequence_query_ext [last_valid_index - (x_max-1) : last_valid_index + 1]
            finestra_mismatch_subject = sequence_sub_ext [last_valid_index - (x_max-1) : last_valid_index + 1]
            #print(finestra_mismatch_query, finestra_mismatch_subject)

            #print(f"Finestra mismatch: {finestra_mismatch_query},{finestra_mismatch_subject}")
            sequences_query, sequences_subject = handle_gaps(finestra_mismatch_query, finestra_mismatch_subject, gap_size,gap_target,x_max)
            #print(sequences_query,sequences_subject)
            if len(sequences_query) > 0:
                for i in range(len(sequences_query)):
                    extension_right_query.append(sequence_query_ext[:last_valid_index - (x_max - 1)] + sequences_query[i])
            #print(f"Estensione query con gap: {extension_right_query}")
            if len(sequences_subject) > 0:
                for i in range(len(sequences_subject)):
                    extension_right_sub.append(sequence_sub_ext[:last_valid_index - (x_max - 1)] + sequences_subject[i])
            #print(f"Estensione subject con gap: {extension_right_subject}")
        else:
            gap_string = "_" * gap_size
            if gap_target == "query":  # Aggiungiamo il gap alla query
                extension_right_query = [sequence_query_ext + gap_string]
                extension_right_sub = [sequence_sub_ext]
            elif gap_target == "subject":  # Aggiungiamo il gap al subject
                extension_right_query = [sequence_query_ext]
                extension_right_sub = [sequence_sub_ext + gap_string]
            
   
    else:
        if last_valid_index >= (x_max - 1):
            extension_right_query = sequence_query_ext[:last_valid_index - (x_max - 1)]
            extension_right_sub = sequence_sub_ext[:last_valid_index - (x_max - 1)]
        else:
            extension_right_query = sequence_query_ext[:last_valid_index + 1]
            extension_right_sub = sequence_sub_ext[:last_valid_index + 1]

        #print(f"Estensione senza gap: {extension_withoutgap_right_query, extension_withoutgap_right_sub}")

    return extension_right_query, extension_right_sub

####ESEMPIO#####
# Sequenze di esempio
sequence_query_ext = "ATGTTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATATCTTGAGTGCAGTTGAGGCAGGC"
sequence_sub_ext =   "TTTTTTAAGTCAGCTGTCAAATGTATGGGCTCAACCCTTATTATCAGTTGAAACTGTAAGACTTGAGTGCGGTGTGGGTACAT"
print(f'lunghezza di partenza = {len(sequence_query_ext)},{len(sequence_sub_ext)}')

# Parametri
#x_max = 6
#perform_gap = False
#gap_size = 2

# Chiamata della funzione
a,b = extend_seed_right(sequence_query_ext, sequence_sub_ext, gap_size=1,perform_gap=True,gap_target='subject',x_max=6)
print("Risultato estensione:", a,b)
print(len(a[0]),len(b[0]))
#score = calculate_score(a[0],b[0])
#print(score)


def extend_seed(df_seeds, query_partenza, list_partenza_subject, x_max):

    contenitore_hsp_query = []
    contenitore_hsp_sub = []
    contenitore_score = []

    hsp_query = ''
    hsp_sub = ''
    score = 0

    query_header = query_partenza[0]
    sequence_query = query_partenza[1]

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

        if len(seeds) > 1:
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

                hsp_query = sequence_query[start_q:end_q]
                hsp_sub = sequence_sub[start_s:end_s]
                score = end_q-start_q

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

                            #print(f'sto valutando questa query: {query_header} con posizione iniziale: {start_q} e finale: {end_q} e la sto confrontando con {col} con posizione iniziale: {start_s} e finale: {end_s} ')
                            
                            sequence_query_ext = sequence_query[end_q:next_start_q]
                            sequence_sub_ext = sequence_sub[end_s:next_start_s]                            
                            extension_right_query,extension_right_sub = extend_seed_right(
                                sequence_query_ext,
                                sequence_sub_ext,
                                gap_size = 0,
                                gap_target = None,
                                perform_gap = False,
                                x_max = 6,
                            )
                            #print(sequence_query[start_q:end_q], sequence_query[next_start_q:next_end_q])
                            #print(sequence_query[start_q:next_end_q])
                            #print(sequence_sub[start_s:end_s], sequence_sub[next_start_s:next_end_s])
                            #print(sequence_sub[start_s:next_end_s])
                            #print(extension_right_query,extension_right_sub)

                            # Aggiorna end_q e end_s
                            new_start_q = end_q
                            new_start_s = end_s
                            new_end_q = end_q + len(extension_right_query[0])
                            new_end_s = end_s + len(extension_right_sub[0])
                            hsp_query += sequence_query[new_start_q:new_end_q]
                            hsp_sub += sequence_sub[new_start_s:new_end_s]
                            #print(hsp_query,hsp_sub)
                            #print(len(hsp_query),len(hsp_sub))
                            score_extension = calculate_score(extension_right_query,extension_right_sub) 
                            #print(score_extension)
                            score += score_extension
                            #print(score)
                            
                        else:
                            # Gestisci i gap
                            if diff_q < diff_s:
                                # Gap nella query
                                gap_size = diff_s - diff_q
                                #print(current_seed)
                                #print(next_seed)
                                #print(sequence_query[:end_q])
                                #print(sequence_sub[:end_s])
                                #print(len(sequence_query[:end_q]))
                                #print(len(sequence_sub[:end_s]))
                                sequence_query_ext = sequence_query[end_q:next_start_q]
                                sequence_sub_ext = sequence_sub[end_s:next_start_s]
                                #print(sequence_query_ext)
                                #print(sequence_sub_ext)
                                #print(len(sequence_query_ext))
                                #print(len(sequence_sub_ext))
                                extension_right_query,extension_right_sub = extend_seed_right(
                                    sequence_query_ext,
                                    sequence_sub_ext,
                                    gap_size,
                                    gap_target='query',
                                    perform_gap=True,
                                    x_max = 6
                                    )
                                #print(len(extension_right_query[0]),len(extension_right_sub[0]))
                                #print(hsp_query,hsp_sub)
                                #print(len(hsp_query),len(hsp_sub))
                                # Aggiorna end_q e end_s
                                hsp_query += extension_right_query[0]
                                hsp_sub += extension_right_sub[0]
                                #print(hsp_query,hsp_sub)
                                #print(len(hsp_query),len(hsp_sub))
                                #print(extension_right_query[0],extension_right_sub[0])
                                score_extension = calculate_score(extension_right_query[0],extension_right_sub[0]) 
                                #print(score_extension)
                                score += score_extension
                                #print(current_seed,next_seed)
                                #print(score)
                                    

                            else:
                                # Gap nel subject
                                gap_size = diff_q - diff_s
                                print(col, current_seed, next_seed)
                                sequence_query_ext = sequence_query[end_q:next_start_q]
                                sequence_sub_ext = sequence_sub[end_s:next_start_s]
                                print(sequence_query_ext)
                                print(sequence_sub_ext)
                                extension_right_query,extension_right_sub = extend_seed_right(
                                    sequence_query_ext,
                                    sequence_sub_ext,
                                    gap_size,
                                    gap_target='subject',
                                    perform_gap=True,
                                    x_max = 6
                                    )
                                print(len(extension_right_query[0]),len(extension_right_sub[0]))
                                print(hsp_query,hsp_sub)
                                print(len(hsp_query),len(hsp_sub))
                                hsp_query += extension_right_query[0]
                                hsp_sub += extension_right_sub[0]
                                score_extension = calculate_score(extension_right_query[0],extension_right_sub[0]) 
                                score += score_extension


                                """
                                # Aggiungi il seed esteso con gap nella query
                                contenitore_seed = []
                                contenitore_score = []
                                for elem in extension_right:
                                    new_end_q = end_q + len(elem)
                                    new_end_s = end_s + len(elem)
                                    extended_seed = [start_q, start_s, new_end_q, new_end_s]
                                    contenitore_seed.append(extended_seed)
                                    new_query = sequence_query[start_q:new_end_q]
                                    new_subject = sequence_sub[start_s:new_end_s]
                                    score_extended_seeds = calculate_score(new_query,new_subject) 
                                    contenitore_score.append(score_extended_seeds)
                                    
                                i = contenitore_score.index(max(contenitore_score))
                                hsp = contenitore_seed[i]
                                """


                else:
                    # Seed finale, potrebbe non avere un successivo
                    hsp_query += sequence_query[start_q:end_q]
                    hsp_sub += sequence_sub[start_s:end_s]
                    score += end_q-start_q
                    contenitore_hsp_query.append(hsp_query)
                    contenitore_hsp_sub.append(hsp_sub)
                    contenitore_score.append(score)
                    break

    return contenitore_hsp_query,contenitore_hsp_sub,contenitore_score
"""
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

#return new_df
"""



#a = extend_seed(seeds1,query_partenza,list_partenza_subject,6)








#p = extend_seed(seeds1,list_partenza_query,list_partenza_subject)
#final_1
#(379,52)
# final_2
# (381,49)



