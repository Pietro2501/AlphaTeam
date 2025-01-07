import argparse
from Query import Query
from Subject import Subject
import ErroriPersonalizzati

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', help='File Query di riferimento', required=True)
    parser.add_argument('-s', '--subject', help='File Subject di riferimento', required=True)
    parser.add_argument('-k','--kmer-length',type=int,default=22,help='Dimensione del kmer da estrarre')
    return parser.parse_args()
'''
query = Query('query.fasta')
query.parse_file()
kmer_query_dict = query.kmer_indexing(22)
kmer_comprev_query_dict = query.kmer_indexing_comp_rev(22)


if kmer_query_dict is None and kmer_comprev_query_dict is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_query_dict,dict) or not isinstance(kmer_comprev_query_dict,dict):
    raise ErroriPersonalizzati.NotADict()

sub = Subject('ref.fa')
sub.parse_file()
kmer_subject_dict = sub.subject_indexing(22)
kmer_comprev_subject_dict = sub.subject_indexing_comp_rev(22)

if kmer_subject_dict is None and kmer_comprev_subject_dict is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_query_dict,dict) or not isinstance(kmer_comprev_query_dict,dict):
    raise ErroriPersonalizzati.NotADict()
    '''

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
                                seed_dict[kmer1]['query'][key1] = [pos1]
                            if key2 not in seed_dict[kmer1]['subject']:
                                seed_dict[kmer1]['subject'][key2] = [pos2]

    for key1,inner_dict in kmer_query_dict.items():
        for kmer1,pos1 in inner_dict.items():
            for key2,sub_dict in kmer_comprev_subject_dict.items():
                for kmer2,pos2 in sub_dict.items():
                    if kmer1 == kmer2:
                        if kmer1 not in seed_dict:
                            seed_dict[kmer1]={'query':{key1:pos1},'subject':{key2:pos2}}
                        else:
                            if key1 not in seed_dict[kmer1]['query']:
                                seed_dict[kmer1]['query'][key1] = [pos1]
                            if key2 not in seed_dict[kmer1]['subject']:
                                seed_dict[kmer1]['subject'][key2] = [pos2]

    return seed_dict

#a = find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)

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





#b = find_comprev_seed(kmer_comprev_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)
#print(a)
#print(b)
#print(len(b))
