import os
import ErroriPersonalizzati
from BioTools import tools,parserFasta



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

    query_partenza e query_partenza_2 : list
        Liste di query generate dalla funzione parse file. Di default impostate a None.

    kmer_query1 : list
        Lista di kmer contenuti nella Query1. Di default impostata a None.

    kmer_query2 : list
        Lista di kmer contenuti nella Query2. Di default impostata a None.

    kmer_query1r : list
        Lista di kmer contenuti nella Query1 complementare revertita. Di default impostata a None.

    kmer_query2r : list
        Lista di kmer contenuti nella Query2 complementare revertita. Di default impostata a None.

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
        self.query_partenza = None
        self.query_partenza_2 = None
        self.kmer_query1 = None
        self.kmer_query2 = None
        self.kmer_query1r = None
        self.kmer_query2r = None
        self.sequence_file = sequence_file
        self.seq_list = None
        self.forward_kmers = None
        self.comp_rev_kmers = None

    def parse_file(self,num_sequences = None):
        """
        Parsa il file sequenza in formato FASTA o nei formati compressi.

        Parameters
        ----------
        num_sequences: int
            Numero di sequenze da processare. Di default impostato a None.
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
        if num_sequences is not None and len(self.seq_list) != num_sequences:
            raise ErroriPersonalizzati.SequenceError(f"Questo BLAST lavora con 2 sequenze ma in realtà ne sono state trovate: {len(self.seq_list)}")

        if num_sequences == 2:
            self.query_partenza = self.seq_list[0]
            self.query_partenza_2 = self.seq_list[1]
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
        if not self.seq_list:
            raise ErroriPersonalizzati.SequenceError('Devi richiamare prima la funzione parse_file()')
        complete_list = []
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
        if not self.seq_list:
            raise ErroriPersonalizzati.SequenceError('Devi richiamare prima la funzione parse_file()')
        complete_list_comp_rev = []
        for element in self.seq_list:
            sequence_comp_rev = tools.fn_comp_rev(element[1])[1]
            kmer_seq_comp_rev = tools.divide_into_kmer(sequence_comp_rev, k)
            complete_list_comp_rev.append(element[0])
            complete_list_comp_rev.append(kmer_seq_comp_rev)
            self.comp_rev_kmers = kmer_seq_comp_rev
        return complete_list_comp_rev

    def slice_for_two_query(self):
        if not self.forward_kmers or not self.comp_rev_kmers:
            raise ErroriPersonalizzati.SequenceError(
                'Devi richiamare prima la funzione di kmer indexing'
            )
        self.kmer_query1 = self.forward_kmers[0:2]
        self.kmer_query2 = self.forward_kmers[2:4]

        self.kmer_query1r = self.comp_rev_kmers[0:2]
        self.kmer_query2r = self.comp_rev_kmers[2:4]
