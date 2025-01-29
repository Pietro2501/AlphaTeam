
import pandas as pd

transizione = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
trasversione = {
    'A': ['C', 'T'],
    'C': ['A', 'G'],
    'G': ['C', 'T'],
    'T': ['A', 'G']
}

def create_query_df(kmer_query_list: list) -> pd.DataFrame:
    """
    La funzione prende una lista di input (kmer_query_list) e costruisce un dataframe Pandas 
    organizzato in modo specifico con intestazioni e indici derivati dalla lista. 
    Parameters
    ----------
     kmer_query_list (list): Lista alternata di intestazioni (index pari) e liste di k-mers (index dispari).
    Return
    ------
     df (pd.DataFrame) : dataframe con relativi valori. 
    """
    headers = []
    kmers = []
    for element in range(0, len(kmer_query_list)):
        if element == 0 or element % 2 == 0:
            headers.append(kmer_query_list[element])
        else:
            kmers.append(kmer_query_list[element])
    df = pd.DataFrame('', index=kmers, columns=headers)
    df = df.reset_index()
    df.columns.values[0] = 'kmer'
    return df
def fill_df_query(df: pd.DataFrame,nome_colonna: str) -> pd.DataFrame:
    """
    La funzione prende in input un dataframe esistente e aggiorna la colonna presa in input, riempiendone le celle
    Parameters
    ----------
     df (pd.DataFrame): dataframe di partenza.
     nome_colonna (str): intestazione della colonna.
    Return
    ------
     df (pd.DataFrame): dataframe aggiornato sulla base della colonna con i suoi valori.
    """
    df[nome_colonna] = range(len(df))
    return df
def create_sub_df(kmer_sub_list: list) -> pd.DataFrame:
    """
    La funzione  costruisce un DataFrame Pandas organizzato in modo specifico con intestazioni e indici derivati dalla lista in input.
    In particolare, gli indici sono unici e disposti in ordine alfabetico.
    Parameters
    ----------
    kmer_sub_list (list): lista alternata di intestazioni delle sequenze e liste di k-mers relativi ad ogni sequenza.
    Return
    ------
    df (pd.DataFrame): dataframe che ha come indici i k-mers unici e in ordine alfabetico e, come nomi delle colonne, le intestazioni delle sequenze subject.
    """
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
def fill_sub_df(df: pd.DataFrame,kmer_subject_list: list) -> pd.DataFrame:
    """
    La funzione prende in input un dataframe pre-esistente e ne riempie le celle ricavando informazioni dalla lista in input.
    In particolare, essa verifica se ogni K-mer è presente nella lista dei k-mers associata ad ogni colonna e inserisce nelle celle 
    l'indice che il k-mer occupa nella lista se esso è presente in quest'ultima.
    Parameters
    ----------
     df (pd.DataFrame): dataframe iniziale
     kmer_subject_list (list): lista alternata di intestazioni delle sequenze e liste di k-mers relativi ad ogni sequenza.
    Return
    ------
     df (pd.DataFrame): dataframe completo di valori.
    """
    for i in range(0,len(kmer_subject_list),2):
        coppia = kmer_subject_list[i:i+2]
        header = coppia[0]
        kmer_list = coppia[1]
        if header in df.columns[1:].tolist():
            for kmer in df['kmer'].tolist():
                if kmer in kmer_list:
                    if kmer_list.count(kmer) == 1:
                        index = kmer_list.index(kmer)
                        df.loc[df['kmer'] == kmer,header] = index
                    elif kmer_list.count(kmer) > 1:
                        raise ValueError()
                else:
                    df.loc[df['kmer'] == kmer,header] = None
    return df
def create_df_with_positions(query_df: pd.DataFrame,subject_df: pd.DataFrame,filename: str) -> pd.DataFrame:
    """
    La funzione combina i due dataframe in input, tramite inner join, per generarne uno nuovo contenente posizioni corrispondenti
    di k-mer tra i due. Salva poi il risultato in un file CSV e restituisce il dataframe finale.
    Parameters
    ----------
     query_df (pd.DataFrame): dataframe con i k-mers realtivi alla sequenza query e relative posizioni.
     subject_df (pd.DataFrame): dataframe con i k-mers realtivi alle sequenze subject e relative posizioni nella lista.
     filename (str): nome del file (senza estensione) in cui il dataframe risultante verrà salvato in formato .csv.
    Return
    ------
     final_df (pd.DataFrame): dataframe risultante che mappa i k-mer con le posizioni della query e dei subject in tuple,
       salvato anche come file CSV.
    """
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
def find_seeds(df: pd.DataFrame,filename:str, kmer_length:int =22) -> pd.DataFrame:
    """
    La funzione identifica regioni di interesse ("seeds") basandosi sulle tuple (query_start, subject_start) e sul valore kmer_length.
    Processa un dataframe contenente tuple rappresentanti posizioni iniziali di k-mer per creare un nuovo dataframe con "seeds"
    che includono le posizioni iniziali e finali. I risultati vengono salvati in un file CSV.
    Parameters
    ----------
     df (pd.DataFrame): dataframe con informazioni realtive alle posizioni dei k-mers in seqeunze query e subject.
     filename (str): nome del file (senza estensione) in cui il dataframe risultante verrà salvato in formato .csv.
     kmer_length (int): lunghezza dei k-mers.
    Return
    ------
     new_df (pd.DataFrame): dataframe con i seeds organizzati, dove ogni seed include posizioni iniziali e finali. 
        I risultati sono salvati anche nel file CSV specificato.
    """
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
                query_end = query_start
                subject_end = subject_start

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


def get_sequence(header: str, list_partenza_subject: list) -> str or None :
    """
    Recupera la sequenza corrispondente a un dato header.
    Parameters
    ----------
     header (str): header del subject.
     list_partenza_subject (list): lista di tuple, ogni tupla composta da header e sequenza associata.
    Return
    ------
    str o None 
    -item[1] (str): se trova corrispondenza
    -None se non trova niente.
    """
    for item in list_partenza_subject:
        if item[0] == header:
            return item[1]
    return None
def calculate_score(sequence_1: str,sequence_2: str) -> int:
    """
    La funzione calcola il punteggio di somiglianza tra due sequenze in base a regole specifiche di corrispondenza e penalità.
    Parameters
    ----------
     sequence_1 (str): prima sequenza da confrontare.
     sequence_2 (str): seconda sequenza da confrontare.
    Return
    ------
     score (int): punteggio di somiglianza.
    """
    score = 0
    stop_calc = min(len(sequence_1),len(sequence_2))
    flag = 0
    for b in range(stop_calc):
        query_char = sequence_1[b]
        sub_char = sequence_2[b]
        
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
    return score
def handle_gaps(finestra_aggiunta_gap_query: str, finestra_aggiunta_gap_sub:str, finestra_mismatch_query: str, finestra_mismatch_sub:str, gap_size: int, gap_target:str, x_max: int) -> tuple[str, str]:
    """
    La funzione gestisce i gap aggiungendoli nella sequenza individuata da gap_target.
    Parameters
    ----------
     finestra_aggiunta_gap_query (str):  sequenza della query dove aggiungere i gap.
     finestra_aggiunta_gap_sub (str): sequenza della subject dove aggiungere i gap.
     finestra_mismatch_query (str): finestra della query che presenta mismatch.
     finestra_mismatch_sub (str): finestra della subject che presenta mismatch.
     gap_size (int): numero massimo di gap consentiti durante l'estensione
     gap_target (str): specifica dove aggiungere i gap. Può essere 'query' o 'subject'.
     x_max (int) : valore massimo negativo del punteggio accettabile per le sequenze.
    Return
    ------
    tuple[str,str]: tupla composta da:
     -finestra_aggiunta_gap_query (str): finestra di sequenza query con aggiunta dei gap, o '*' se i gap non migliorano il punteggio.
     -finestra_aggiunta_gap_sub (str): finestra di sequenza subject con aggiunta dei gap, o '*' se i gap non migliorano il punteggio.
    """
    max_gap = 3
    sequence_query_list = []
    sequence_sub_list = []
    score_list = []
    gap_numbers = min(gap_size,max_gap) # viene calcolato in modo da capire quanti gap inserire, per limitarne l'espansione

    if gap_target == 'query':
        for a in range(gap_numbers):
            sequence_finestra_gap = '_' * (a + 1) + finestra_mismatch_query
            gap_score = calculate_score(sequence_finestra_gap, finestra_mismatch_sub)
            if gap_score >= -x_max:
                finestra_aggiunta_gap_query = '_' * (a + 1) + finestra_aggiunta_gap_query
                gap_score = calculate_score(finestra_aggiunta_gap_query, finestra_aggiunta_gap_sub)
                sequence_query_list.append(finestra_aggiunta_gap_query)
                score_list.append(gap_score)
        if len(score_list) > 0 and len(sequence_query_list) > 0:
            maximum = max(score_list)
            i = score_list.index(maximum)
            finestra_aggiunta_gap_query = sequence_query_list[i]
        else:
            finestra_aggiunta_gap_query = '*'

    if gap_target == 'subject':
        for a in range(gap_numbers):
            sequence_finestra_gap = '_' * (a + 1) + finestra_mismatch_sub
            gap_score = calculate_score(sequence_finestra_gap, finestra_mismatch_query)
            if gap_score >= -x_max:
                finestra_aggiunta_gap_sub  = '_' * (a + 1) + finestra_aggiunta_gap_sub
                gap_score = calculate_score(finestra_aggiunta_gap_query , finestra_aggiunta_gap_sub)
                sequence_sub_list.append(finestra_aggiunta_gap_sub )
                score_list.append(gap_score)
        if len(score_list) > 0 and len(sequence_sub_list) > 0:
            maximum = max(score_list)
            i = score_list.index(maximum)
            finestra_aggiunta_gap_sub = sequence_sub_list[i]
        else:
            finestra_aggiunta_gap_sub = '*'
    
    return finestra_aggiunta_gap_query,finestra_aggiunta_gap_sub
def find_mismatch_window(sequence_query_ext: str, sequence_sub_ext:str, x_max: int) -> tuple[int,int]:
    """
    La funzione confronta le basi nelle stesse posizioni delle sequenze in input assegnando un punteggio positivo ad ogni match.
    Ogni mismatch consecutivo incrementa un contatore. Se il numero di mismatch consecutivi raggiunge il valore soglia,
    il ciclo si interrompe.
    Parameters
    ----------
     sequence_query_ext (str): sequenza estesa da query da confrontare.
     sequence_sub_ext (str): sequenza estesa da subject da confrontare.
     x_max (int): numero massimo di mismatch consecutivi consentiti prima di interrompere il confronto.
    Return
    ------
    tuple[int,int]: tupla composta da:
     -last_valid_index (int): indice dell'ultima posizione valida prima dei mismatch consecutivi
     -mismatch_consecutivi (int): numero di mismatch consecutivi rilevati.
    """
    mismatch_consecutivi = 0
    last_valid_index = 0

    for a in range(len(sequence_query_ext)):
        if a >= len(sequence_sub_ext):
            break

        if sequence_query_ext[a] == sequence_sub_ext[a]:
            mismatch_consecutivi = 0
        else:
            mismatch_consecutivi += 1

            if mismatch_consecutivi == x_max:
                last_valid_index = a
                break

    return last_valid_index, mismatch_consecutivi
def extend_seed_right(sequence_query_ext:str, sequence_sub_ext:str, gap_size:int,perform_gap:bool,gap_target:str, x_max:int) -> tuple[str,str]: 
    """
    Estende l'allineamento della query e delle sequenze subject da destra, inserisce i gap quando è necessario.
    Parameters
    ----------
    sequence_query_ext (str): sequenza query da estendere.
    sequence_sub_ext (str): sequenza subject da estendere.
    gap_size (int): numero massimo di gaps consentiti durante l'estensione.  
    perform_gap (bool): se True, i gaps sono inseriti per ottimizzare l'allineamento quando i mismatches superano la soglia consentita.
    gap_target (str):specifica dove introdurre i gaps: “query” aggiunge i gaps alla sequenza query e “subject” li aggiunge alla sequenza subject. 
    x_max (int): threshold per i mismatches consecutivi necessari per attivare l'inserimento dei gaps o fermare l'estensione.  
    Return
    ------
    tuple[str,str]: tupla composta da:
    -extension_right_query (str): sequenza query estesa. 
    -extension_right_sub (str): sequenze subject estese.
    """
    gap_used = 0
    if perform_gap:

        while gap_used < gap_size:
            last_valid_index, mismatch_consecutivi = find_mismatch_window(sequence_query_ext,sequence_sub_ext,x_max)

            if mismatch_consecutivi == x_max:
                finestra_mismatch_query = sequence_query_ext [last_valid_index - (x_max-1) : last_valid_index + 1]
                finestra_mismatch_subject = sequence_sub_ext [last_valid_index - (x_max-1) : last_valid_index + 1]

                finestra_aggiunta_gap_query = sequence_query_ext [last_valid_index - (x_max-1) : ]
                finestra_aggiunta_gap_sub = sequence_sub_ext [last_valid_index - (x_max-1) : ]

                gap_rimasti = gap_size-gap_used
                sequence_query, sequence_subject = handle_gaps(finestra_aggiunta_gap_query, finestra_aggiunta_gap_sub, finestra_mismatch_query, finestra_mismatch_subject, gap_rimasti,gap_target,x_max)

                if sequence_query == '*' or sequence_subject == '*':
                    extension_right_query = sequence_query_ext[:last_valid_index - (x_max - 1)]
                    extension_right_sub = sequence_sub_ext[:last_valid_index - (x_max - 1)]
                else:
                    sequence_query_ext = sequence_query_ext[:last_valid_index - (x_max - 1)] + sequence_query


                    sequence_sub_ext = sequence_sub_ext[:last_valid_index - (x_max - 1)] + sequence_subject

                    extension_right_query = sequence_query_ext
                    extension_right_sub = sequence_sub_ext

            else:
                gap_string = "_" * (gap_size-gap_used)
                if gap_target == "query":
                    extension_right_query = sequence_query_ext + gap_string
                    extension_right_sub = sequence_sub_ext
                elif gap_target == "subject":
                    extension_right_query = sequence_query_ext
                    extension_right_sub = sequence_sub_ext + gap_string

            gap_aggiunti = extension_right_query.count('_') + extension_right_sub.count('_')
            gap_used += gap_aggiunti

    else:
        last_valid_index, mismatch_consecutivi = find_mismatch_window(sequence_query_ext,sequence_sub_ext,x_max)
        if mismatch_consecutivi == x_max:
            extension_right_query = sequence_query_ext[:last_valid_index - (x_max - 1)]
            extension_right_sub = sequence_sub_ext[:last_valid_index - (x_max - 1)]
        elif mismatch_consecutivi < x_max:
            extension_right_query = sequence_query_ext
            extension_right_sub = sequence_sub_ext

    return extension_right_query, extension_right_sub
def extend_seed(df_seeds:pd.DataFrame, query_partenza: tuple, list_partenza_subject: list, x_max:int) -> tuple[list,list,list,list]:  
    """
    Estende gli allineamenti dei seed iterando su di essi in un DataFrame fra la sequenza query e le molteplici sequenze subject per produrre degli high-scoring segment pairs (HSPs).  
    Gestisce i mismatches, i gaps, e gli allineamenti continui laddove possibile.
    Parameters
    ----------
    df_seeds (pd.DataFrame): dataframe contenente i seeds per ogni colonna delle sequenze dei subject. 
    query_partenza (tuple): tupla contenente l'header della query e la sua sequenza.
    list_partenza_subject (list): lista di tuple, in cui ogni tupla contiene l'header della sequenza subject e la sequenza corrispondente. 
    x_max (int): numero massimo di mismatches consecutivi consentito prima dell'arresto dell'estensione dell'allineamento. 
    Return
    ------
    tuple[list,list,list,list]: tupla composta da:
    -contenitore_hsp_query (list): lista degli allineamenti della query. 
    -contenitore_hsp_sub (list): lista degli allineamenti delle subject estese. 
    -contenitore_score (list): lista degli scores per ogni allineamento esteso, in cui lo score si basa sulla lunghezza e sui matches.  
    -contenitore_ref (list): lista degli header dei subject del corrispettivo allineamento esteso.  
    """
    contenitore_hsp_query = []
    contenitore_hsp_sub = []
    contenitore_score = []
    contenitore_ref = []

    query_header = query_partenza[0]
    sequence_query = query_partenza[1]

    extended_seeds_dict = {}

    for col in df_seeds.columns:
        seeds = df_seeds[col].dropna().tolist()

        extended_seeds = []

        sequence_sub = get_sequence(col, list_partenza_subject)
        if sequence_sub is None:
            extended_seeds_dict[col] = extended_seeds
            continue

        list_inner_hsp_q = []
        list_inner_hsp_s = []
        list_inner_score = []

        if len(seeds) > 0:
            for i, current_seed in enumerate(seeds):
                if not current_seed:
                    continue

                start_q, start_s, end_q, end_s = map(int, current_seed)

                hsp_query = sequence_query[start_q:end_q]
                hsp_sub = sequence_sub[start_s:end_s]
                score = end_q - start_q

                if i < len(seeds)-1:
                    next_seed = seeds[i + 1]
                    if next_seed:
                        next_start_q, next_start_s, next_end_q, next_end_s = next_seed

                        try:
                            next_start_q = int(next_start_q)
                            next_start_s = int(next_start_s)
                            next_end_q = int(next_end_q)
                            next_end_s = int(next_end_s)
                        except ValueError:
                            print(f"Errore: I valori di next_start_q, next_start_s, next_end_q, next_end_s devono essere interi. Seed: {next_seed}")
                            continue

                        diff_q = next_start_q - end_q
                        diff_s = next_start_s - end_s

                        if diff_q == diff_s:
                            sequence_query_ext = sequence_query[end_q:next_start_q]
                            sequence_sub_ext = sequence_sub[end_s:next_start_s]

                            extension_right_query, extension_right_sub = extend_seed_right(
                                sequence_query_ext,
                                sequence_sub_ext,
                                gap_size=0,
                                gap_target=None,
                                perform_gap=False,
                                x_max=6,
                            )

                        else:
                            if diff_q < diff_s:
                                gap_size = diff_s - diff_q

                                sequence_query_ext = sequence_query[end_q:next_start_q]
                                sequence_sub_ext = sequence_sub[end_s:next_start_s]

                                extension_right_query,extension_right_sub = extend_seed_right(
                                    sequence_query_ext,
                                    sequence_sub_ext,
                                    gap_size,
                                    gap_target='query',
                                    perform_gap=True,
                                    x_max = 6
                                    )

                            else:
                                gap_size = diff_q - diff_s

                                sequence_query_ext = sequence_query[end_q:next_start_q]
                                sequence_sub_ext = sequence_sub[end_s:next_start_s]

                                extension_right_query,extension_right_sub = extend_seed_right(
                                    sequence_query_ext,
                                    sequence_sub_ext,
                                    gap_size,
                                    gap_target='subject',
                                    perform_gap=True,
                                    x_max = 6
                                    )

                        hsp_query += extension_right_query
                        hsp_sub += extension_right_sub

                        score_extension = calculate_score(extension_right_query,extension_right_sub)
                        score += score_extension
                        list_inner_hsp_q.append(hsp_query)
                        list_inner_hsp_s.append(hsp_sub)
                        list_inner_score.append(score)

                else:
                    list_inner_hsp_q.append(hsp_query)
                    list_inner_hsp_s.append(hsp_sub)
                    list_inner_score.append(score)

        hsp_query_completo = ''.join(list_inner_hsp_q)
        hsp_sub_completo = ''.join(list_inner_hsp_s)
        score_completo = sum(list_inner_score)

        contenitore_hsp_query.append(hsp_query_completo)
        contenitore_hsp_sub.append(hsp_sub_completo)
        contenitore_score.append(score_completo)
        contenitore_ref.append(col)
    return contenitore_hsp_query, contenitore_hsp_sub, contenitore_score, contenitore_ref


def create_results_table(cont_hsp_query: list, cont_hsp_sub: list, cont_hsp_score: list, col: list,filename: str) -> pd.DataFrame: 
   """
   Crea una result table dagli HSP (High-scoring Segment Pair) data e li salva in un file .csv. 
   La tabella include query sequence, subject sequences, headers e gli scores, ordinati in modo decrescente. 
   Parameters
   ----------
   cont_hsp_query (list): lista delle sequenze query derivate dall'analisi dell'HSP.
   cont_hsp_sub (list): lista delle sequenze subject derivate dall'analisi dell'HSP.
   cont_hsp_score (list): lista di score corrispondente agli HSP.
   col (list): elenco di headers o identificatori per i risultati HSP.
   filename (str):nome del file (senza estensione) nel quale il result-table verrà salvato.
   Return
   ------
   results_df (pd.DataFrame): salvato come file .csv.Un pandas DataFrame contenente
     i risultati degli HSP con le colonne 'hsp_query','hsp_sub','header','hsp_score'.
   """ 
   data = {'hsp_query':cont_hsp_query,
                'hsp_sub':cont_hsp_sub,
                'header': col,
                'hsp_score':cont_hsp_score,
                }
   results_df = pd.DataFrame(data)
   results_df = results_df.sort_values(['hsp_score'], ascending=False)
   results_df.to_csv(filename + '.csv')
   return results_df
def metrics_for_blast(best_alignment_df:pd.DataFrame,query_partenza: tuple,list_partenza_subject:list) -> tuple[list,list,list]:
    """
    La funzione calcola le metriche degli allineamenti del BLAST, includendo query_coverage, E-value e la percentuale di identità. 
    Parameters
    ----------
    best_alignment_df (pd.Dataframe): dataframe contenente gli allineamenti migliori.  
    query_partenza (tuple): tupla in cui il secondo elemento è una sequenza query usata per le analisi del BLAST.  
    list_partenza_subject (list): lista di sequenze subject che recupera la lunghezza dei subject per il calcolo dell' E-value. 
    Return
    ------
    tuple[list,list,list]: tupla composta da:
    -query_cov_list (list): lista della copertura in percentuale della query per ogni allineamento.
    -e_value_list (list): lista degli E-values calcolato per ogni allineamento.
    -identity_list (list): lista dei valori in percentuale dell'identità per ogni allineamento.
    """
    #coverage = lung allineamento / lung query * 100
    #E-value
    #E_value = len(query) * len(sbj_seqs) * 2 ** (-score)
    # Identità = quanti sono i match rispetto all'allineamento
    query_cov_list = []
    e_value_list = []
    identity_list = []

    query_col = best_alignment_df.columns[0]
    sub_col = best_alignment_df.columns[1]
    score_col = best_alignment_df.columns[3]
    header_col = best_alignment_df.columns[2]
    for q,s in zip(best_alignment_df[query_col].tolist(),best_alignment_df[sub_col].tolist()):
        match = 0
        if len(q) == len(s):
            query_coverage = len(q) / len(query_partenza[1]) * 100
            query_cov_list.append(query_coverage)

            for base in range(len(q)):
                query_char = q[base]
                sub_char = s[base]
                if query_char == sub_char:
                    match += 1
                else:
                    continue

            max_identity = round(match/len(q)*100,2)
            identity_list.append(max_identity)

    for score,col in zip(best_alignment_df[score_col].tolist(),best_alignment_df[header_col].tolist()):
        len_sub = len(get_sequence(col,list_partenza_subject))
        len_query = len(query_partenza[1])
        e_value = len_query * len_sub * 2 ** (-score)
        e_value_list.append(e_value)

    return query_cov_list,e_value_list,identity_list
def blast_result_df(best_alignment_query:pd.DataFrame,query_partenza:tuple,list_partenza_subject:list,filename:str) -> pd.DataFrame:
    """
    La funzione migliora un dataframe di allineamento BLAST con le metriche calcolate e le salva in un DataFrame come file .csv.
    Parameters
    ----------
    best_alignment_query (pd.DataFrame): dataframe contenente i migliori allineamenti.
    query_partenza (tuple): tupla in cui il secondo elemento della query è la sequenza query usata per le analisi del BLAST. 
    list_partenza_subject (list): lista di sequenze subject che ritorna la lunghezza delle subject per il calcolo dell’E-value.
    filename (str): nome del file nel quale il dataframe verrà salvato.
    Return
    ------
    best_alignment_query (pd.DataFrame): dataframe salvato in un file .csv con lo specifico filename.
    """
    coverage_cont, e_value_cont, identity_cont = metrics_for_blast(best_alignment_query,query_partenza,list_partenza_subject)
    best_alignment_query = best_alignment_query.drop(['hsp_query','hsp_sub'],axis = 1)
    best_alignment_query.insert(loc=2, column='query_coverage', value=coverage_cont)
    best_alignment_query.insert(loc=3, column='E_value', value=e_value_cont)
    best_alignment_query.insert(loc=4, column='max_identity', value=identity_cont)

    best_alignment_query.to_csv(filename + '.csv')

    return best_alignment_query
def print_alignment(best_alignment_df:pd.DataFrame,output_file :str,result_df:pd.DataFrame) -> None:
    """
    La funzione rappresenta i risultati degli allineamenti in un file, includendo l'allineamento di sequenza, le metriche e le annotazioni.
    Parameters
    ----------
    best_alignment_df (pd.DataFrame): dataframe contenente i migliori allineamenti in cui nelle colonne inseriamo le query e le sequenze delle subject.
    output_file (str): percorso del file di output in cui verranno scritti i risultati dell'allineamento formattato.
    result_df (pd.DataFrame): dataframe contenente le metriche calcolate per ogni allineamento
    Return
    ------
    None
        Scrive i dettagli dell'allineamento nello specifico output finale.
    """
    query_hsp = best_alignment_df.columns[0]
    sub_hsp = best_alignment_df.columns[1]

    with open(output_file, 'w') as f:
        for (index,row1), (index,row2) in zip(best_alignment_df.iterrows(),result_df.iterrows()):
            q = row1[query_hsp]
            s= row1[sub_hsp]
            sub_header = row2.iloc[0]
            score = row2.iloc[1]
            query_coverage = row2.iloc[2]
            e_value = row2.iloc[3]
            max_identity = row2.iloc[4]
            f.write(
                f"Subject: {sub_header}, Score: {score}, Coverage: {query_coverage:.2f}%, E-value: {e_value:.2e}, Max_identity:{max_identity:.2f}%\n")

            middle_line = []
            length = max(len(q), len(s))
            for base in range(length):
                query_char = q[base] if base < len(q) else '-'
                sub_char = s[base] if base < len(s) else '-'
                if query_char == '_' or sub_char == '_' or query_char == '-' or sub_char == '-':
                    middle_line.append(' ')
                elif query_char == sub_char:
                    middle_line.append('|')
                else:
                    middle_line.append('.')

            middle_str = ''.join(middle_line)

            f.write(f"Query:   {q}\n")
            f.write(f"         {middle_str}\n")
            f.write(f"Subject: {s}\n")
            f.write('\n')









