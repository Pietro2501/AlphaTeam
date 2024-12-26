from BioTools import tools


# OopCompanion:suppressRename

class Query:
    def __init__(self, query_file):
        self.query_file = query_file

    def parse_file(self):
        if self.query_file.endswith('.gz') or self.query_file.endswith('.gzip'):
            tools.extract_info(self.query_file, 'query.txt')
        else:
            print("Il file è già nell'estensione .txt!")
            return self.query_file

    def kmer_indexing(self, k: int) -> tuple:
        with open(self.query_file, 'r') as seq:
            seq = seq.readlines()
            seq = ''.join(seq)
        kmer_set = set([seq[i:i + k] for i in range(len(seq) - k + 1)])
        self.kmer_set = kmer_set
        return kmer_set

    def kmer_indexing_comp_rev(self, k: int) -> tuple:
        with open(self.query_file, 'r') as seq:
            seq = seq.readlines()
            seq = ''.join(seq)
            seq_comp_rev = tools.fn_comp_rev(seq)[1]
        kmer_set_comp_rev = set([seq_comp_rev[i:i + k] for i in range(len(seq_comp_rev) - k + 1)])
        self.kmer_set_comp_rev = kmer_set_comp_rev
        return kmer_set_comp_rev


query = Query('query.txt')
query.parse_file()
print("Stampo i kmer della query")
print(query.kmer_indexing(11))
print("Stampo i kmer della query del complementare revertito")
print(query.kmer_indexing_comp_rev(11))