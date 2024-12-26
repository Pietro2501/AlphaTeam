from BioTools import tools
from BioTools import parserFasta


# OopCompanion:suppressRename
class Subject:
    def __init__(self,subject_file):
        self.subject_file = subject_file

    def parse_file(self):
        if self.subject_file.endswith('.gz') or self.subject_file.endswith('.gzip'):
            extracted_file = tools.extract_info(self.subject_file, 'subject')
            self.subject_file = extracted_file
        else:
            print("Il file è già nell'estensione .fasta!")
            self.subject_file = self.subject_file

    def subject_indexing(self):
        diz = parserFasta.parse_fasta(self.subject_file)
        kmer_diz = {}
        k = 22
        for header,sequence in diz.items():
            kmers = tools.kmer_counting(sequence,k)
            kmer_diz[header] = {}
            for i,kmer in enumerate(kmers):
                kmer_diz[header][kmer] = i*k
        return kmer_diz


sub = Subject('esempioFasta.fasta')
sub.parse_file()
print(sub.subject_indexing())


