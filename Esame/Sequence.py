import os
import ErroriPersonalizzati
from settings import nucleotides_scoring_matrix
import tools, tools1
import parserFasta


# OopCompanion:suppressRename



class Sequence:
    def __init__(self, sequence_file,scoring_matrix = None):
        if not os.path.isfile(sequence_file):
           raise ErroriPersonalizzati.FileNotFoundError(f"File non trovato: {sequence_file}")
        self.sequence_file = sequence_file
        if scoring_matrix is None:
            self.scoring_matrix = nucleotides_scoring_matrix
        else:
            self.scoring_matrix = scoring_matrix
        self.diz = None
        self.forward_kmers = None
        self.comp_rev_kmers = None

    def parse_file(self):
        self.list = list(parserFasta.parse_fasta(self.sequence_file).items())

        if self.sequence_file.endswith('.gz') or self.sequence_file.endswith('.gzip'):
            try:
                tools.extract_info(self.sequence_file, 'query.txt')
            except Exception as e:
                raise ErroriPersonalizzati.SequenceError(f"Errore durante l'estrazione del file compresso: {e}")
        elif not self.sequence_file.endswith('.fasta') and not self.sequence_file.endswith('.fa'):
            raise ErroriPersonalizzati.FileTypeError()
        return self.list

    def kmer_indexing(self, k: int) -> list:
        try:
            #diz = parserFasta.parse_fasta(self.query_file)
            complete_list = []
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for element in self.list:
            #print(element[1])
            kmer_seq = tools.divide_into_kmer(element[1], k)
            complete_list.append(element[0])
            complete_list.append(kmer_seq)
        self.forward_kmers = complete_list
        return complete_list

    def kmer_indexing_comp_rev(self, k: int) -> list:
        try:
            #diz = parserFasta.parse_fasta(self.query_file)
            complete_list_comp_rev = []
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for element in self.list:
            #print(type(abc))
            sequence_comp_rev = tools1.fn_comp_rev(element[1])[1]
            kmer_seq_comp_rev = tools.divide_into_kmer(sequence_comp_rev,k)
            complete_list_comp_rev.append(element[0])
            complete_list_comp_rev.append(kmer_seq_comp_rev)
            self.comp_rev_kmers = kmer_seq_comp_rev
        return complete_list_comp_rev






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

#query = Sequence('query.fasta')
#query.parse_file()
#print(query.kmer_indexing(22))
#print(len(query.kmer_indexing(22)))
#print(query.kmer_indexing_comp_rev(22))