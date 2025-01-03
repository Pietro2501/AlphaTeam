#import sys
#print(sys.path)
#sys.path.append('c:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam')
#import os
#os.chdir ('c:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam')
from BioTools import parserFasta
from BioTools import tools
import ErroriPersonalizzati
import os


def divide_into_kmer(sequence,k):
    """
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
    """
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

# OopCompanion:suppressRename
class Subject:
    def __init__(self,subject_file):
        if not os.path.isfile(subject_file):
            raise ErroriPersonalizzati.FileNotFoundError()
        self.subject_file = subject_file

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
        if self.subject_file.endswith('.gz') or self.subject_file.endswith('.gzip'):
            try:
                extracted_file = tools.extract_info(self.subject_file, 'subject')
                self.subject_file = extracted_file
            except Exception as e:
                raise ErroriPersonalizzati.Error(f"File non trovato:{self.subject_file}")
        #elif not self.subject_file.endswith('.fasta'):
            #raise ErroriPersonalizzati.FileTypeError()
        else:
            print("Il file è già nell'estensione .fasta!")
            self.subject_file = self.subject_file

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

    def subject_indexing(self):
        """
        Effettua il parsing del file FASTA e divide ogni sequenza in kmer di lunghezza 22.

        Parametri:
        Lavora sulla variabile di istanza presente nel self

        Return:
            complete_dict: dict
            Dizionario in cui le chiavi sono gli header del Fasta e i valori
            sono dizionari di kmer di lunghezza 22, con le posizioni in cui ciascun k-mer appare.
        """
        try:
            diz = parserFasta.parse_fasta(self.subject_file)
            complete_dict = diz
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for header, sequence in complete_dict.items():
            complete_dict[header] = divide_into_kmer(sequence,22)
        return complete_dict

    def subject_indexing_comp_rev(self):
        """
                Effettua il parsing del file FASTA e divide ogni sequenza in kmer di lunghezza 22 e
                applica il complementare revertito della sequenza, prima di suddividerla.

                Parametri:
                Lavora sulla variabile di istanza presente nel self

                Return:
                    complete_dict: dict
                    Dizionario in cui le chiavi sono gli header del Fasta e i valori
                    sono dizionari di kmer di lunghezza 22 , con le posizioni in cui ciascun k-mer appare.
                """
        try:
            diz = parserFasta.parse_fasta(self.subject_file)
            complete_dict = diz
        except Exception as e:
            raise ErroriPersonalizzati.FastaParsingError()
        for header, sequence in complete_dict.items():
            #CHIEDERE A MELANIA
            complete_dict[header] = divide_into_kmer(tools.fn_comp_rev(sequence)[1],22)
        return complete_dict

sub = Subject('/Users/pietrodispaldro/Documents/GitHub/AlphaTeam/esempioFASTA.fa')
sub.parse_file()
print(sub.subject_indexing())
print(sub.subject_indexing_comp_rev())



