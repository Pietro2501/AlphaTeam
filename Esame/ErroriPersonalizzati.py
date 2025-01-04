
class Error(Exception):
    pass

class FileTypeError(Error):
    def __init__(self, message="Il file deve avere estensione .fa, .fasta, .gz, o .gzip"):
        super().__init__(message)

class FileNotFoundError(Error):
    def __init__(self, message="File non trovato o percorso errato"):
        super().__init__(message)

class KmerError(Error):
    def __init__(self, message="Il valore di k deve essere un intero positivo"):
        super().__init__(message)

class KmerTooLong(Error):
    def __init__(self,message="Il valore di k è superiore alla lunghezza della sequenza"):
        super().__init__(message)

class FastaParsingError(Error):
    def __init__(self, message="Errore nel parsing del file FASTA"):
        super().__init__(message)

class QueryError(Error):
    def __init__(self, message="Errore durante l'estrazione del file"):
        super().__init__(message)