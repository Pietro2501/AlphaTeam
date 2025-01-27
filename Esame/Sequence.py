import os
import ErroriPersonalizzati
import tools, tools1
import parserFasta


# OopCompanion:suppressRename

class Sequence:
    """
    Una classe che gestisce e processa i file di sequenza e 
    indicizza i k-mer per le sequenze forward e le loro complementari revertite. 
    Attributes
    ----------
    sequence_file: str
        Il percorso del file della sequenza di input (.fasta, .fa, or compressed .gz).

    seq_list: list
        Parsa le sequenze come liste di tuple. Di default è None. 

    forward_kmers: list
        Sequenze indicizzate di k-mer del filamento forward. Di default è None.

    comp_rev_kmers: list
        Sequenze indicizzate di k-mer del filamento complementare revertito. Di default è None.
    Methods
    -------
        parse_file():
            Ritorna un file parsato. 
        kmer_indexing(k: int) -> list:
            Genera i kmer indicizati per le sequenze forward.
        kmer_indexing_comp_rev(k: int) -> list:
            Genera i kmer indicizati per le sequenze complementari reverite.
    """

    def __init__(self, sequence_file):
        """
        Costruisce tutti gli attributi necessari per l'oggetto Sequence. 
        Parameters
        ----------
        sequence_file: str
            Il percorso (path) del file della sequenza di input (.fasta, .fa, or compressed .gz).
        seq_list: list
            Parsa le sequenze come lista di tuple (ID, sequence). Di default è None.
        forward_kmers: list 
            Sequenze indicizzate di k-mer del filamento forward. Di default è None.
        comp_rev_kmers: list
            Sequenze indicizzate di k-mer del filamento complementare revertito. Di default è None.
        """
        if not os.path.isfile(sequence_file):
           raise ErroriPersonalizzati.FileNotFoundError(f"File non trovato: {sequence_file}")
        self.sequence_file = sequence_file
        self.seq_list = None
        self.forward_kmers = None
        self.comp_rev_kmers = None

    def parse_file(self):
        """
        Parsa il file sequenza in formato FASTA o nei formati compressi.
        Return
        ------
        seq_list: list
            Una lista di sequenze parsate[(ID, sequence)] 
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
        Genera indici k-mer per le sequenze forward.
        Args
        ----
        k: int
            la lunghezza dei k-mers.   
        Return
        ------
        complete_list: list 
            Una lista contenente gli ID e le loro sequenze k-mer corrispondenti.
        """
        try:
            complete_list = []
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for element in self.seq_list:
            kmer_seq = tools.divide_into_kmer(element[1], k)
            complete_list.append(element[0])
            complete_list.append(kmer_seq)
        self.forward_kmers = complete_list
        return complete_list

    def kmer_indexing_comp_rev(self, k: int) -> list:
        """
        Genera indici k-mer per le sequenze complementari revertite. 
        Args
        ----
        k: int
            La lunghezza dei k-mers. 
        Return
        ------
        complete_list_comp_rev: list
            Una lista contenente gli ID e le loro sequenze k-mer complementari revertite corrispondenti.
        """
        try:
            complete_list_comp_rev = []
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for element in self.seq_list:
            sequence_comp_rev = tools1.fn_comp_rev(element[1])[1]
            kmer_seq_comp_rev = tools.divide_into_kmer(sequence_comp_rev,k)
            complete_list_comp_rev.append(element[0])
            complete_list_comp_rev.append(kmer_seq_comp_rev)
            self.comp_rev_kmers = kmer_seq_comp_rev
        return complete_list_comp_rev
