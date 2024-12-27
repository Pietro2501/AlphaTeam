def divide_into_kmer(sequence,k):
    kmer_diz = {}
    for i in range(0, len(sequence)-k+1,k):
        kmer = sequence[i:i+k]
        kmer_diz.setdefault(kmer, [])
        kmer_diz[kmer].append(i)
    return kmer_diz

#seq = 'ACGTGGACCTTTGAGACGACGTGG'
#a = divide_into_kmer(seq,3)
#print(a)

def subject_indexing(diz):
    complete_dict = diz
    for header, sequence in complete_dict.items():
        complete_dict[header] = divide_into_kmer(sequence,3)
    return complete_dict

diz = {'Homo Sapiens': 'ACGTGGACCTTTGAGACGACGTGG', 'Mus Musculus': 'ACGTTT'}
b = subject_indexing(diz)
print(b)


#
#len(kmer_diz)