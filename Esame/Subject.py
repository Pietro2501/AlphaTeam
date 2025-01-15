import sys
#print(sys.path)
#sys.path.append('C:\Users\\verak\\OneDrive - Università degli Studi di Bari\\Documenti\\GitHub\\AlphaTeam\\Esame')
#os.chdir ('C:\Users\verak\OneDrive - Università degli Studi di Bari\\Documenti\\GitHub\\AlphaTeam\\Esame')
#from BioTools import parserFasta
#from BioTools import tools
import ErroriPersonalizzati
import os
#print(os.path.exists("C:/Users/verak/OneDrive - Università degli Studi di Bari/Documenti/GitHub/AlphaTeam/Esame"))
import tools
import tools1 as tools1
import parserFasta

"""
def divide_into_kmer(sequence,k):

    Divide la sequenza in kmer di lunghezza k, restituendo un dizionario
    dove la chiave sono i kmer e il valore è una lista di posizioni.

    Verranno sollevati errori personalizzati se:
    - k è minore o uguale a zero
    - k è maggiore della lunghezza della sequenza

    Parametri:
    sequence: str
        Sequenza di caratteri
    k: int
        Dimensione fissa del kmer da estrarre a partire dalla sequenza

    Return:
    kmer_diz: dict
        Un dizionario in cui le chiavi sono i kmer e i valori sono
        le liste di posizioni.

    if k <= 0:
        raise ErroriPersonalizzati.KmerError()
    if k > len(sequence):
        raise ErroriPersonalizzati.KmerTooLong()
    kmer_diz = {}
    for i in range(0, len(sequence)-k+1,k):
        kmer = sequence[i:i+k]
        kmer_diz.setdefault(kmer, [])
        kmer_diz[kmer].append(i)
    return kmer_diz
"""

# OopCompanion:suppressRename
class Subject:
    def __init__(self,subject_file):
        if not os.path.isfile(subject_file):
            raise ErroriPersonalizzati.FileNotFoundError()
        self.subject_file = subject_file
        self.sub_diz = None
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
                If the file is not in a valid FASTA format (e.g., wrong extension).
            """
        self.sub_diz = parserFasta.parse_fasta(self.subject_file)
        if self.subject_file.endswith('.gz') or self.subject_file.endswith('.gzip'):
            try:
                extracted_file = tools1.extract_info(self.subject_file,'Subject')
                self.subject_file = extracted_file
            except Exception as e:
                raise ErroriPersonalizzati.Error(f"File non trovato:{self.subject_file}")
        elif not self.subject_file.endswith('.fa') or self.subject_file.endswith('.fasta'):
            raise ErroriPersonalizzati.FileTypeError()
        else:
            print("Il file è già nell'estensione .fasta/.fa!")
        
        return self.sub_diz


    #def subject_indexing(self):
    #    diz = parserFasta.parse_fasta(self.subject_file)
    #    kmer_diz = {}
    #    k = 22
    #    for header,sequence in diz.items():
    #        kmers = tools.kmer_counting(sequence,k)
    #        kmer_diz[header] = {}
    #        for i,kmer in enumerate(kmers):
    #            kmer_diz[header][kmer] = i*k
    #    return kmer_diz

    def subject_indexing(self,k):
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
            #diz = parserFasta.parse_fasta(self.subject_file)
            complete_dict_sub = {}
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for header, sequence in self.sub_diz.items():
            complete_dict_sub[header] = tools.divide_into_kmer(sequence,k)
        self.forward_kmers = complete_dict_sub
        return complete_dict_sub

    def subject_indexing_comp_rev(self,k):
        """"
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
            #diz = parserFasta.parse_fasta(self.subject_file)
            complete_dict_comp_rev = {}
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for header, sequence in self.sub_diz.items():
            #CHIEDERE A MELANIA
            complete_dict_comp_rev[header] = tools.divide_into_kmer(tools1.fn_comp_rev(sequence)[1],k)
            self.comp_rev_kmers = complete_dict_comp_rev
        return complete_dict_comp_rev

#sub = Subject('C:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam\\Esame\\ref.fa')
#sub.parse_file()
#print(sub.subject_indexing(22))
#print(sub.subject_indexing_comp_rev(22))