import os
import ErroriPersonalizzati
from settings import nucleotides_scoring_matrix
import tools, tools1
import parserFasta


# OopCompanion:suppressRename

class Sequence:
    """
    A class for handling and processing sequence files and 
    generate k-mer indices for forward sequences and their complementary reverse. 
    Attributes
    ----------
    sequence_file: str
        The path to the input sequence file (.fasta, .fa, or compressed .gz).

    seq_list: list
        Parsed sequences as a list of tuples. Default is None.

    forward_kmers: list
        Kmer indexed sequences from the forward strand. Default is None.

    comp_rev_kmers: list
        Kmer indexed sequences from the complementary reverse strand. Default is None.
    Methods
    -------
        parse_file():
            returns a parsed file.
        kmer_indexing(k: int) -> list:
            Generates kmer indices for forward sequences.
        kmer_indexing_comp_rev(k: int) -> list:
            Generates k-mer indices for complementary reverse sequences
    """

    def __init__(self, sequence_file):
        """
        Constructs all the necessary attributes for the Sequence object.
        Parameters
        ----------
        sequence_file: str
            The path to the input sequence file (.fasta, .fa, or compressed .gz).
        seq_list: list
            Parsed sequences as a list of tuples (ID, sequence). Default is None.
        forward_kmers: list 
            K-mer indexed sequences from the forward strand. Default is None.
        comp_rev_kmers: list
            K-mer indexed sequences from the complementary reverse strand. Default is None.
        Raise
        -----
            FileNotFoundError: If the specified file does not exist.
        """
        if not os.path.isfile(sequence_file):
           raise ErroriPersonalizzati.FileNotFoundError(f"File non trovato: {sequence_file}")
        self.sequence_file = sequence_file
        self.seq_list = None
        self.forward_kmers = None
        self.comp_rev_kmers = None

    def parse_file(self):
        """
        Parses the sequence file in FASTA and compressed files.
        Return
        ------
        seq_list: list
            A list of parsed sequences in the format [(ID, sequence)]       
        Raise
        -----
        SequenceError: If an error occurs while extracting a compressed file.
        FileTypeError: If the file is not a valid FASTA format (.fasta, .fa).
        """
        if self.sequence_file.endswith('.gz') or self.sequence_file.endswith('.gzip'):
            try:
                tools.extract_info(self.sequence_file, 'query.txt')
            except Exception as e:
                raise ErroriPersonalizzati.SequenceError(f"Errore durante l'estrazione del file compresso: {e}")
        elif not self.sequence_file.endswith('.fasta') and not self.sequence_file.endswith('.fa'):
            raise ErroriPersonalizzati.FileTypeError()
        self.seq_list = list(parserFasta.parse_fasta(self.sequence_file).items())
        return self.seq_list

    def kmer_indexing(self, k: int) -> list:
        """
        Generates k-mer indices for forward sequences.
        Args
        ----
        k: int
            The length of the k-mers.   
        Return
        ------
        complete_list: list 
            A list containing IDs and their corresponding k-mer sequences.
        Raise
        -----
        FastaParsingError: If an error occurs while parsing the FASTA file.
        """
        try:
            complete_list = []
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for element in self.list:
            kmer_seq = tools.divide_into_kmer(element[1], k)
            complete_list.append(element[0])
            complete_list.append(kmer_seq)
        self.forward_kmers = complete_list
        return complete_list

    def kmer_indexing_comp_rev(self, k: int) -> list:
        """
        Generates k-mer indices for complementary reverse sequences.
        Args
        ----
        k: int
            The length of the k-mers.
        Return
        ------
        complete_list_comp_rev: list
            A list containing IDs and their corresponding complementary reverse k-mer sequences.
        Raise
        -----
        FastaParsingError: If an error occurs while parsing the FASTA file.
        """
        try:
            complete_list_comp_rev = []
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for element in self.list:
            sequence_comp_rev = tools1.fn_comp_rev(element[1])[1]
            kmer_seq_comp_rev = tools.divide_into_kmer(sequence_comp_rev,k)
            complete_list_comp_rev.append(element[0])
            complete_list_comp_rev.append(kmer_seq_comp_rev)
            self.comp_rev_kmers = kmer_seq_comp_rev
        return complete_list_comp_rev
