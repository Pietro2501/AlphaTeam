
import Blast, Sequence,ErroriPersonalizzati,os
from BioTools import tools

output_folder = 'blast_results'
os.makedirs(output_folder, exist_ok=True)


query = Sequence.Sequence('query.fasta')
list_partenza_query = query.parse_file(2)
kmer_query_list = query.kmer_indexing(22)
kmer_comprev_query_list = query.kmer_indexing_comp_rev(22)
if kmer_query_list is None and kmer_comprev_query_list is None:
    raise ErroriPersonalizzati.EmptyList()
if not isinstance(kmer_query_list, list) or not isinstance(kmer_comprev_query_list, list):
    raise ErroriPersonalizzati.NotAList()


query.slice_for_two_query()

sub = Sequence.Sequence('ref.fa')
list_partenza_subject = sub.parse_file()
kmer_subject_list = sub.kmer_indexing(22)
kmer_comprev_subject_list = sub.kmer_indexing_comp_rev(22)

########## CREAZIONE QUERY DF #############
dataf1 = Blast.create_query_df(query.kmer_query1)
dataf2 = Blast.create_query_df(query.kmer_query2)
datafr1 = Blast.create_query_df(query.kmer_query1r)
datafr2 = Blast.create_query_df(query.kmer_query1r)

############## FASE DI FILL DEI QUERY DF ################
dataf1_fill = Blast.fill_df_query(dataf1, "b6635d67cb594473ddba9f8cfba5d13d")
dataf2_fill = Blast.fill_df_query(dataf2, "4516aa60a483dd8c7bbc57098c45f1a5")
datafr1_fill = Blast.fill_df_query(datafr1, "b6635d67cb594473ddba9f8cfba5d13d")
datafr2_fill = Blast.fill_df_query(datafr2, "4516aa60a483dd8c7bbc57098c45f1a5")

########## CREAZIONE SUB DF #############
datasub = Blast.create_sub_df(kmer_subject_list)

############## CARICAMENTO DF SUB FILLATO
data = tools.load_table('filled_sub_df.csv')

################ TABELLE PER RICERCA SEED #########################
df_final_1 = Blast.create_df_with_positions(dataf1_fill, data, os.path.join(output_folder,'final1_df'))
df_final_2 = Blast.create_df_with_positions(dataf2_fill, data, os.path.join('final2_df'))

seeds1 = Blast.find_seeds(df_final_1, os.path.join(output_folder,'query1seeds'))
seeds2 = Blast.find_seeds(df_final_2, os.path.join(output_folder,'query2seeds'))

hsp_query, hsp_sub, score, col = Blast.extend_seed(seeds1, query.query_partenza, list_partenza_subject, 6)
hsp_query2, hsp_sub2, score2, col2 = Blast.extend_seed(seeds2, query.query_partenza_2, list_partenza_subject, 6)


all_alignments_query1 = Blast.create_results_table(hsp_query, hsp_sub, score, col, os.path.join(output_folder,'all_alignments_query1'))

all_alignments_query2 = Blast.create_results_table(hsp_query2, hsp_sub2, score2, col2, os.path.join(output_folder,'all_alignments_query2'))

best_alignment_query_1 = all_alignments_query1.iloc[:10].copy()
best_alignment_query_2 = all_alignments_query2.iloc[:10].copy()

worst_alignment_query_1 = all_alignments_query1.iloc[-10:].copy()
worst_alignment_query_2 = all_alignments_query2.iloc[-10:].copy()

metrics_query_1 = Blast.metrics_for_blast(best_alignment_query_1, query.query_partenza, list_partenza_subject)
metrics_query_2 = Blast.metrics_for_blast(best_alignment_query_2, query.query_partenza_2, list_partenza_subject)

coverage_cont, e_value_cont, identity_cont = metrics_query_2


metrics_worst_query_1 = Blast.metrics_for_blast(worst_alignment_query_1, query.query_partenza, list_partenza_subject)
metrics_worst_query_2 = Blast.metrics_for_blast(worst_alignment_query_2, query.query_partenza_2, list_partenza_subject)

blast_result1 = Blast.blast_result_df(best_alignment_query_1, query.query_partenza, list_partenza_subject, os.path.join(output_folder,'blast_result1'))
blast_result2 = Blast.blast_result_df(best_alignment_query_2, query.query_partenza_2, list_partenza_subject, os.path.join(output_folder,'blast_result2'))

blast_worst_result1 = Blast.blast_result_df(worst_alignment_query_1, query.query_partenza, list_partenza_subject,
                                      os.path.join(output_folder,'blast_worst_result1'))
blast_worst_result2 = Blast.blast_result_df(worst_alignment_query_2, query.query_partenza_2, list_partenza_subject,
                                      os.path.join(output_folder,'blast_worst_result2'))

Blast.print_alignment(best_alignment_query_1, os.path.join(output_folder,'allineamento_top10_query1.txt'), blast_result1)
Blast.print_alignment(best_alignment_query_2, os.path.join(output_folder,'allineamento_top10_query2.txt'), blast_result2)

Blast.print_alignment(worst_alignment_query_1, os.path.join(output_folder,'allineamento_flop10_query1.txt'), blast_worst_result1)
Blast.print_alignment(worst_alignment_query_2, os.path.join(output_folder,'allineamento_flop10_query2.txt'), blast_worst_result2)










