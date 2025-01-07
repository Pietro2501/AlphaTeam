from Query import Query
from Subject import Subject
import ErroriPersonalizzati

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

b = find_comprev_seed(kmer_comprev_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)
print(a)
print(b)
#print(len(b))
