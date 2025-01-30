import traceback
import Blast,Sequence,ErroriPersonalizzati,os,argparse,time
import sys
from BioTools import tools

output_folder = 'blast_results'
os.makedirs(output_folder, exist_ok=True)

def main():
    start_time = time.time()
    try:
        print("\033[1;34m🔄 Avvio del processo BLAST...\033[0m")
        time.sleep(1)
        print("\033[1;34m📌 Parsing degli argomenti...\033[0m")
        time.sleep(1)
        parser = argparse.ArgumentParser(description="BLAST script")
        parser.add_argument("-q", "--query", required=True, help="Path della Query")
        parser.add_argument("-s", "--subject", required=True, help="Path della Subject")
        parser.add_argument("-k", "--kmer_size", required=True,type=int, help="Lunghezza k-mer")
        args = parser.parse_args()

        print("\033[1;34m📂 Caricamento della query...\033[0m")
        time.sleep(0.5)
        query = Sequence.Sequence(args.query)
        query.parse_file(2)
        print("\033[1;32m✅ Query caricata con successo!\033[0m")
        time.sleep(1)

        print("\033[1;34m🔍 Indicizzazione dei k-mer della query...\033[0m")
        time.sleep(0.5)
        kmer_query_list = query.kmer_indexing(22)
        kmer_comprev_query_list = query.kmer_indexing_comp_rev(22)
        if kmer_query_list is None and kmer_comprev_query_list is None:
            raise ErroriPersonalizzati.EmptyList()
        if not isinstance(kmer_query_list, list) or not isinstance(kmer_comprev_query_list, list):
            raise ErroriPersonalizzati.NotAList()
        query.slice_for_two_query()
        print("\033[1;32m✅ Indicizzazione k-mer completata!\033[0m")
        time.sleep(1)

        print("\033[1;34m📂 Caricamento del subject...\033[0m")
        time.sleep(0.5)
        sub = Sequence.Sequence(args.subject)
        list_partenza_subject = sub.parse_file()
        kmer_subject_list = sub.kmer_indexing(22)
        kmer_comprev_subject_list = sub.kmer_indexing_comp_rev(22)
        print("\033[1;32m✅ Subject caricato con successo!\033[0m")
        time.sleep(1)

        ########## CREAZIONE QUERY DF #############
        print("\033[1;34m📊 Creazione e Fill DataFrame della query...\033[0m")
        time.sleep(0.5)
        dataf1 = Blast.create_query_df(query.kmer_query1,"b6635d67cb594473ddba9f8cfba5d13d")
        dataf2 = Blast.create_query_df(query.kmer_query2,"4516aa60a483dd8c7bbc57098c45f1a5")
        datafr1 = Blast.create_query_df(query.kmer_query1r,"b6635d67cb594473ddba9f8cfba5d13d")
        datafr2 = Blast.create_query_df(query.kmer_query1r,"4516aa60a483dd8c7bbc57098c45f1a5")
        print("\033[1;32m✅ DataFrame query creati e riempiti con successo!\033[0m")
        time.sleep(1)

        ########## CREAZIONE SUB DF #############
        print("\033[1;34m📊 Creazione DataFrame del subject...\033[0m")
        time.sleep(1)
        datasub = Blast.create_sub_df(kmer_subject_list)
        print("\033[1;32m✅ DataFrame del subject creati con successo\033[0m")
        time.sleep(1)

        ############## CARICAMENTO DF SUB FILLATO
        data = tools.load_table('filled_sub_df.csv')

        ################ TABELLE PER RICERCA SEED #########################
        df_final_1 = Blast.create_df_with_positions(dataf1, data, os.path.join(output_folder, 'final1_df'))
        df_final_2 = Blast.create_df_with_positions(dataf2, data, os.path.join(output_folder,'final2_df'))

        print("\033[1;34m🔍 Ricerca dei seed nei DataFrame...\033[0m")
        time.sleep(0.5)
        seeds1 = Blast.find_seeds(df_final_1, os.path.join(output_folder, 'query1seeds'))
        seeds2 = Blast.find_seeds(df_final_2, os.path.join(output_folder, 'query2seeds'))
        print("\033[1;32m✅ Ricerca seed completata!\033[0m")
        time.sleep(1)

        print("\033[1;34m🔄 Estensione dei seed...\033[0m")
        time.sleep(0.5)
        hsp_query, hsp_sub, score, col = Blast.extend_seed(seeds1, query.query_partenza, list_partenza_subject, 6)
        hsp_query2, hsp_sub2, score2, col2 = Blast.extend_seed(seeds2, query.query_partenza_2, list_partenza_subject, 6)
        print("\033[1;32m✅ Estensione dei seed completata!\033[0m")
        time.sleep(1)

        print("\033[1;34m📑 Generazione e salvataggio dei risultati...\033[0m")
        time.sleep(0.5)
        all_alignments_query1 = Blast.create_results_table(hsp_query, hsp_sub, score, col, os.path.join(output_folder, 'all_alignments_query1'))
        all_alignments_query2 = Blast.create_results_table(hsp_query2, hsp_sub2, score2, col2, os.path.join(output_folder, 'all_alignments_query2'))

        best_alignment_query_1 = all_alignments_query1.iloc[:10].copy()
        best_alignment_query_2 = all_alignments_query2.iloc[:10].copy()

        worst_alignment_query_1 = all_alignments_query1.iloc[-10:].copy()
        worst_alignment_query_2 = all_alignments_query2.iloc[-10:].copy()

        metrics_query_1 = Blast.metrics_for_blast(best_alignment_query_1, query.query_partenza, list_partenza_subject)
        metrics_query_2 = Blast.metrics_for_blast(best_alignment_query_2, query.query_partenza_2, list_partenza_subject)

        metrics_worst_query_1 = Blast.metrics_for_blast(worst_alignment_query_1, query.query_partenza, list_partenza_subject)
        metrics_worst_query_2 = Blast.metrics_for_blast(worst_alignment_query_2, query.query_partenza_2, list_partenza_subject)

        blast_result1 = Blast.blast_result_df(best_alignment_query_1, query.query_partenza, list_partenza_subject, os.path.join(output_folder, 'blast_result1'))
        blast_result2 = Blast.blast_result_df(best_alignment_query_2, query.query_partenza_2, list_partenza_subject, os.path.join(output_folder, 'blast_result2'))

        blast_worst_result1 = Blast.blast_result_df(worst_alignment_query_1, query.query_partenza, list_partenza_subject,
                                                    os.path.join(output_folder,'blast_worst_result1'))
        blast_worst_result2 = Blast.blast_result_df(worst_alignment_query_2, query.query_partenza_2, list_partenza_subject,
                                                    os.path.join(output_folder,'blast_worst_result2'))


        Blast.print_alignment(best_alignment_query_1, os.path.join(output_folder, 'allineamento_top10_query1.txt'), blast_result1)
        Blast.print_alignment(best_alignment_query_2, os.path.join(output_folder, 'allineamento_top10_query2.txt'), blast_result2)

        Blast.print_alignment(worst_alignment_query_1, os.path.join(output_folder, 'allineamento_flop10_query1.txt'), blast_worst_result1)
        Blast.print_alignment(worst_alignment_query_2, os.path.join(output_folder, 'allineamento_flop10_query2.txt'), blast_worst_result2)
        print("\033[1;32m✅ Risultati salvati con successo!\033[0m")
        tempo_passato = time.time() -start_time

        print("\033[1;32m🎉 Esecuzione completata con successo!\033[0m")
        print(f'\033[1;32mIl main è stato eseguito con successo in {tempo_passato: .2f} secondi \033[0m')
        print(f'\033[1;32mTrovi i risultati nella cartella: {output_folder}\033[0m')
        return True
    except Exception as e:
        print("\033[1;31m❌ Il programma ha incontrato un errore!\033[0m")
        print(f"\033[1;31mErrore: {e}\033[0m")
        traceback.print_exc()
        return False

if __name__ == '__main__':
    main()











