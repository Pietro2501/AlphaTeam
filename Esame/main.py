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

            if current_value is not None and isinstance(current_value,tuple):

                if prev_value is None and next_value is not None:
                    query_start, subject_start = current_value
                    column_seeds.append(((query_start, subject_start), None))

                elif prev_value is not None and next_value is None:
                    query_end, subject_end = current_value
                    query_end += kmer_length
                    subject_end += kmer_length

                    if column_seeds and column_seeds[-1][1] is None:
                        column_seeds[-1] = (column_seeds[-1][0], (query_end, subject_end))


                elif prev_value is None and next_value is None:
                    query_start, subject_start = current_value
                    query_end = query_start + kmer_length
                    subject_end = subject_start + kmer_length
                    column_seeds.append(((query_start, subject_start), (query_end, subject_end)))


        seed_results.append(column_seeds)


    max_rows = max(len(seeds) for seeds in seed_results)
    new_df_data = {df.columns[i]: seed_results[i] + [None] * (max_rows - len(seed_results[i])) for i in range(len(df.columns))}
    new_df = pd.DataFrame(new_df_data)

    new_df.index = range(1, len(new_df) + 1)
    new_df.to_csv(filename +'.csv')
    return new_df
seeds1 = find_seeds(df_final_1,'query1seeds')
seeds2 = find_seeds(df_final_2,'query2seeds')



#final_1
#(379,52)
# final_2
# (381,49)



