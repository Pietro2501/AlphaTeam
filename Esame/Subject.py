#import sys
#print(sys.path)
#sys.path.append('c:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam')
#import os
#os.chdir ('c:\\Users\\Melania\\Documents\\GitHub\\AlphaTeam')
from BioTools import tools
from BioTools import parserFasta

def divide_into_kmer(sequence,k):
    kmer_diz = {}
    for i in range(0, len(sequence)-k+1,k):
        kmer = sequence[i:i+k]
        kmer_diz.setdefault(kmer, [])
        kmer_diz[kmer].append(i)
    return kmer_diz

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
        diz = parserFasta.parse_fasta(self.subject_file)
        complete_dict = diz
        for header, sequence in complete_dict.items():
            complete_dict[header] = divide_into_kmer(sequence,22)
        return complete_dict

    def subject_indexing_comp_rev(self):
        diz = parserFasta.parse_fasta(self.subject_file)
        complete_dict = diz
        for header, sequence in complete_dict.items():
            complete_dict[header] = divide_into_kmer(tools.fn_comp_rev(sequence)[1],22)
        return complete_dict

sub = Subject('esempioFASTA.fa')
sub.parse_file()
print(sub.subject_indexing())
print(sub.subject_indexing_comp_rev())



