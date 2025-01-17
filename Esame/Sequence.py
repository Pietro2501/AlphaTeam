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
        """
            Parses a query file, processes it, and handles compressed files.

            This method reads a query file in FASTA format, processes it using a FASTA parser,
            and handles cases where the file is compressed. If the file is not in the expected format,
            it raises custom exceptions.

            Returns:
            dict
                A dictionary where keys are sequence identifiers and values are the corresponding sequences
                from the parsed FASTA file.

            Raises:
            ErroriPersonalizzati.QueryError
                If an error occurs during the extraction of a compressed file.

            ErroriPersonalizzati.FileTypeError
                If the file is not in a valid FASTA format
            """
        self.diz = parserFasta.parse_fasta(self.sequence_file)

        if self.sequence_file.endswith('.gz') or self.sequence_file.endswith('.gzip'):
            try:
                tools.extract_info(self.sequence_file, 'query.txt')
            except Exception as e:
                raise ErroriPersonalizzati.SequenceError(f"Errore durante l'estrazione del file compresso: {e}")
        elif not self.sequence_file.endswith('.fasta') and not self.sequence_file.endswith('.fa'):
            raise ErroriPersonalizzati.FileTypeError()
        return self.diz

    def kmer_indexing(self, k: int) -> dict:
        """
            Generates k-mer indexes for all sequences in the parsed FASTA file.

            This method divides each sequence in the parsed FASTA dictionary into k-mers of specified length
            and returns a dictionary mapping sequence headers to their respective k-mers.

            Parameters:
            k : int
                The length of the k-mers to generate.

            Returns:
            dict
                A dictionary where keys are sequence headers and values are lists of k-mers generated
                from the corresponding sequences.

            Raises:
            ErroriPersonalizzati.FastaParsingError
                If an error occurs while parsing the FASTA file.
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
            Generates k-mer indexes for the reverse complement of all sequences in the parsed FASTA file.

            This method computes the reverse complement of each sequence in the parsed FASTA dictionary,
            divides it into k-mers of specified length, and returns a dictionary mapping sequence headers
            to their respective k-mers.

            Parameters:
            k : int
                The length of the k-mers to generate.

            Returns:
            dict
                A dictionary where keys are sequence headers and values are lists of k-mers generated
                from the reverse complement of the corresponding sequences.

            Raises:
            ErroriPersonalizzati.FastaParsingError
                If an error occurs while parsing the FASTA file.
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

query = Sequence('query.fasta')
query.parse_file()
print(query.kmer_indexing(22))