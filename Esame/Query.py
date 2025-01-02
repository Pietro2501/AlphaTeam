import os
import ErroriPersonalizzati
from BioTools import tools
from settings import nucleotides_scoring_matrix
import tools




# OopCompanion:suppressRename

class Query:
    def __init__(self, query_file,scoring_matrix = None):
        if not os.path.isfile(query_file):
           raise ErroriPersonalizzati.FileNotFoundError(f"File non trovato: {query_file}")
        self.query_file = query_file
        if scoring_matrix is None:
            self.scoring_matrix = nucleotides_scoring_matrix
        else:
            self.scoring_matrix = scoring_matrix

    def parse_file(self):
        """
                Analizza il file. Se questo è compresso (formato .gz o .gzip)
                ,ne tenta l'estrazione. Se il file non è in formato .fasta,
                solleva l'eccezione corrispondente.

                Parametri:
                Lavora sulla variabile di istanza presente nel self

                Return:
                Nessuno.
        """
        if self.query_file.endswith('.gz') or self.query_file.endswith('.gzip'):
            try:
                tools.extract_info(self.query_file, 'query.txt')
            except Exception as e:
                raise ErroriPersonalizzati.QueryError(f"Errore durante l'estrazione del file compresso: {e}")
        elif not self.query_file.endswith('.txt'):
            raise ErroriPersonalizzati.FileTypeError()
        else:
            print("Il file è già nell'estensione .txt!")
            return self.query_file

    def kmer_indexing(self, k: int) -> set:
        """
                 Divide la query in kmer di lunghezza 22.

                Parametri:
                Lavora sulla variabile di istanza presente nel self
                    k: int
                    Dimensione fissa del kmer da estrarre

                Return:
                    kmer_set: set
                    Set contenente i kmer estratti dalla sequenza query
                """
        if k <= 0:
            raise ErroriPersonalizzati.KmerError()
        with open(self.query_file, 'r') as seq:
            seq = seq.readlines()
            seq = ''.join(seq)
            if k > len(seq):
                raise ErroriPersonalizzati.KmerTooLong()
        kmer_set = set([seq[i:i + k] for i in range(len(seq) - k + 1)])
        self.kmer_set = kmer_set
        return kmer_set

    def kmer_indexing_comp_rev(self, k: int) -> set:
        """
                         Divide il complementare revertito della query
                         in kmer di lunghezza fissa 22.

                        Parametri:
                        Lavora sulla variabile di istanza presente nel self
                            k: int
                            Dimensione fissa del kmer da estrarre

                        Return:
                            kmer_set_com_rev: set
                            Set contenente i kmer estratti dalla sequenza query,
                            lavorando sul complementare revertito
                        """
        if k <= 0:
            raise ErroriPersonalizzati.KmerError()
        with open(self.query_file, 'r') as seq:
            seq = seq.readlines()
            seq = ''.join(seq)
            if k > len(seq):
                raise ErroriPersonalizzati.KmerTooLong()
            seq_comp_rev = tools.fn_comp_rev(seq)[1]
        kmer_set_comp_rev = set([seq_comp_rev[i:i + k] for i in range(len(seq_comp_rev) - k + 1)])
        self.kmer_set_comp_rev = kmer_set_comp_rev
        return kmer_set_comp_rev

    def generate_word_diz(self,k:int,threshold = 0,max_words = None) -> dict:
        #if not self.kmer_set:
           #self.kmer_indexing(k)

        words_diz = {}
        for kmer in self.kmer_set:
            words_list = tools.generate_words(
                kmer = kmer,
                scoring_matrix= self.scoring_matrix,
                threshold=threshold,
                max_words=max_words
            )
            words_diz[kmer] = words_list

        self.words_diz = words_diz
        return words_diz




query = Query('query.txt')
query.parse_file()
print("Stampo i kmer della query")
print(query.kmer_indexing(11))
print("Stampo i kmer della query del complementare revertito")
#print(query.kmer_indexing_comp_rev(11))

print(query.generate_word_diz(22,20,10))