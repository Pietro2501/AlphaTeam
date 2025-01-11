import os
import ErroriPersonalizzati
from settings import nucleotides_scoring_matrix
import tools, tools1
import parserFasta




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
        self.diz = None
        self.forward_kmers = None
        self.comp_rev_kmers = None

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
        self.diz = parserFasta.parse_fasta(self.query_file)

        if self.query_file.endswith('.gz') or self.query_file.endswith('.gzip'):
            try:
                tools.extract_info(self.query_file, 'query.txt')
            except Exception as e:
                raise ErroriPersonalizzati.QueryError(f"Errore durante l'estrazione del file compresso: {e}")
        elif not self.query_file.endswith('.fasta'):
            raise ErroriPersonalizzati.FileTypeError()
        return self.diz

    def kmer_indexing(self, k: int) -> dict:
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
        try:
            #diz = parserFasta.parse_fasta(self.query_file)
            complete_dict = {}
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for header, sequence in self.diz.items():
            complete_dict[header] = tools.divide_into_kmer(sequence,k)
        self.forward_kmers = complete_dict
        return complete_dict


    def kmer_indexing_comp_rev(self, k: int) -> dict:
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
        try:
            #diz = parserFasta.parse_fasta(self.query_file)
            complete_dict_comp_rev = {}
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for header,sequence in self.diz.items():
            #print(type(abc))
            sequence_comp_rev = tools1.fn_comp_rev(sequence)[1]
            complete_dict_comp_rev[header] = tools.divide_into_kmer(sequence_comp_rev,k)
            self.comp_rev_kmers = complete_dict_comp_rev
        return complete_dict_comp_rev

    """
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
    """




#query = Query('query.fasta')
#print(f"prova{query.parse_file()}")
#query.kmer_indexing(22)
#print(query.complete_dict)
#print("Stampo i kmer della query")
#print(query.kmer_indexing(22))
#print("Stampo i kmer della query del complementare revertito")
#print(query.kmer_indexing_comp_rev(22))

#print(f"ciao{query.prova()}")

#print(query.generate_word_diz(22,20,10))