


# OopCompanion:suppressRename

class FileTypeError(Exception):
    def __init__(self, message="Il file deve avere estensione .fa, .fasta, .gz, o .gzip"):
        super().__init__(message)

class FileNotFoundError(Exception):
    def __init__(self, message="File non trovato o percorso errato"):
        super().__init__(message)

class KmerError(Exception):
    def __init__(self, message="Il valore di k deve essere un intero positivo"):
        super().__init__(message)

class KmerTooLong(Exception):
    def __init__(self,message="Il valore di k è superiore alla lunghezza della sequenza"):
        super().__init__(message)

class FastaParsingError(Exception):
    def __init__(self, message="Errore nel parsing del file FASTA"):
        super().__init__(message)

class SequenceError(Exception):
    def __init__(self, message="Errore durante l'estrazione del file"):
        super().__init__(message)

class EmptyDict(Exception):
    def __init__(self, message="Il dizionario è vuoto"):
        super().__init__(message)

class NotADict(Exception):
    def __init__(self, message="Questo non è un dizionario"):
        super().__init__(message)