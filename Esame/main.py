import argparse
from Query import Query
from Subject import Subject
import ErroriPersonalizzati
"""
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', help='File Query di riferimento', required=True)
    parser.add_argument('-s', '--subject', help='File Subject di riferimento', required=True)
    parser.add_argument('-k','--kmer-length',type=int,default=22,help='Dimensione del kmer da estrarre')
    return parser.parse_args()

"""
query = Query('Esame\\query.fasta')
query.parse_file()
kmer_query_dict = query.kmer_indexing(22)
kmer_comprev_query_dict = query.kmer_indexing_comp_rev(22)


if kmer_query_dict is None and kmer_comprev_query_dict is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_query_dict,dict) or not isinstance(kmer_comprev_query_dict,dict):
    raise ErroriPersonalizzati.NotADict()

sub = Subject('Esame\\ref.fa')
sub.parse_file()
kmer_subject_dict = sub.subject_indexing(22)
kmer_comprev_subject_dict = sub.subject_indexing_comp_rev(22)

if kmer_subject_dict is None and kmer_comprev_subject_dict is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_query_dict,dict) or not isinstance(kmer_comprev_query_dict,dict):
    raise ErroriPersonalizzati.NotADict()

def find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)->dict:

    seed_dict = {}
    for key1,inner_dict in kmer_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_dict:
                            seed_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_dict[kmer1]['query']:
                                seed_dict[kmer1]['query'][key1] = pos1
                            else:
                                for p in pos1:
                                        if p not in seed_dict[kmer1]['query'][key1]:
                                            seed_dict[kmer1]['query'][key1].append(p)

                            if key2 not in seed_dict[kmer1]['subject']:
                                seed_dict[kmer1]['subject'][key2] = pos2
                            else:
                                for p in pos2:
                                    if p not in seed_dict[kmer1]['subject'][key2]:
                                        seed_dict[kmer1]['subject'][key2].append(p)

    for key1,inner_dict in kmer_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_comprev_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_dict:
                            seed_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_dict[kmer1]['query']:
                                seed_dict[kmer1]['query'][key1] = pos1
                            else:
                                for p in pos1:
                                        if p not in seed_dict[kmer1]['query'][key1]:
                                            seed_dict[kmer1]['query'][key1].append(p)

                            if key2 not in seed_dict[kmer1]['subject']:
                                seed_dict[kmer1]['subject'][key2] = pos2
                            else:
                                for p in pos2:
                                    if p not in seed_dict[kmer1]['subject'][key2]:
                                        seed_dict[kmer1]['subject'][key2].append(p)
                        

    return seed_dict

a = find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)

def find_comprev_seed(kmer_comprev_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)->dict:
    seed_comprev_dict = {}
    for key1,inner_dict in kmer_comprev_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_comprev_dict:
                            seed_comprev_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_comprev_dict[kmer1]['query']:
                                seed_comprev_dict[kmer1]['query'][key1] = [pos1]
                            if key2 not in seed_comprev_dict[kmer1]['subject']:
                                seed_comprev_dict[kmer1]['subject'][key2] = [pos2]

    for key1,inner_dict in kmer_comprev_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_comprev_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_comprev_dict:
                            seed_comprev_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_comprev_dict[kmer1]['query']:
                                seed_comprev_dict[kmer1]['query'][key1] = [pos1]
                            if key2 not in seed_comprev_dict[kmer1]['subject']:
                                seed_comprev_dict[kmer1]['subject'][key2] = [pos2]

    return seed_comprev_dict

#def extend_seed(seed_dict,query,subject):
print(a)
#print(len(a))

"""
s = 38
x = 6
schema = []
for kmer, inner_dict in a.items():
    for que_sub, sub_dict in inner_dict.items():
        if que_sub == 'query':
            for j in sub_dict.items():
                query = j
        if que_sub == 'subject':
            for i in sub_dict.items():
                subject = i
                schema.append(kmer)
                schema.append(query)
                schema.append(subject)

#print(schema)

lista_utile=[]
for i,j in kmer_query_dict['b6635d67cb594473ddba9f8cfba5d13d'].items():
    if j == [176]:
        lista_utile.append(i)
seq_utile = ''.join(lista_utile) 
print(seq_utile)

lista_confronto=[]
for x,y in kmer_subject_dict['MJ030-2-barcode67-umi101484bins-ubs-3'].items():
    if y == [506]:
        lista_confronto.append(x)
seq_confronto = ''.join(lista_confronto)
print(seq_confronto)

transizione = {'A':'G','G':'A','C':'T','T':'C'}
trasversione = {'A':'C','A':'T','C':'A','C':'G','G':'C','G':'T','T':'A','T':'G'}

score=22
x=0

for a in range(0,len(seq_utile)):
    if seq_utile[a] == seq_confronto[a]:
        score += 1
        print(a)
        print(seq_utile[a])
        print(seq_confronto[a])
        print(score)
    else:
        chiave = seq_utile[a]
        if transizione[chiave] == seq_confronto[a]:
            x+=1
        if trasversione[chiave] == seq_confronto[a]:
            score -= 1
            x+=1
        print(a)
        print(seq_utile[a])
        print(seq_confronto[a])
        print(score)


print(score)

        

"""




#return

"""
def main():
    args = parse_args()

    query = Query(args.query)
    query.parse_file()

    subject = Subject(args.subject)
    subject.parse_file()

    kmer_query_dict = query.kmer_indexing(args.kmer_length)
    kmer_comprev_query_dict = query.kmer_indexing_comp_rev(args.kmer_length)

    if kmer_query_dict is None and kmer_comprev_query_dict is None:
        raise ErroriPersonalizzati.EmptyDict()
    if not isinstance(kmer_query_dict, dict) or not isinstance(kmer_comprev_query_dict, dict):
        raise ErroriPersonalizzati.NotADict()

    kmer_subject_dict = subject.subject_indexing(args.kmer_length)
    kmer_comprev_subject_dict = subject.subject_indexing_comp_rev(args.kmer_length)

    if kmer_subject_dict is None and kmer_comprev_subject_dict is None:
        raise ErroriPersonalizzati.EmptyDict()
    if not isinstance(kmer_subject_dict, dict) or not isinstance(kmer_comprev_subject_dict, dict):
        raise ErroriPersonalizzati.NotADict()

    a = find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)
    b = find_comprev_seed(kmer_comprev_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)

    print(f"Seed del forward: \033[1m{a}\033[0m")
    print()
    print(f"Seed del revertito complementare: \033[1m{b}\033[0m")

if __name__ == "__main__":
    main()



"""
#b = find_comprev_seed(kmer_comprev_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)

#print(len(a))
#print(b)
#print(len(b))
