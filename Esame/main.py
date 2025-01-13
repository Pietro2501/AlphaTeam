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
query = Query('/Users/pietrodispaldro/Documents/GitHub/AlphaTeam/Esame/query.fasta')
diz_partenza_query =query.parse_file()
kmer_query_dict = query.kmer_indexing(22)
kmer_comprev_query_dict = query.kmer_indexing_comp_rev(22)


if kmer_query_dict is None and kmer_comprev_query_dict is None:
    raise ErroriPersonalizzati.EmptyDict()
if not isinstance(kmer_query_dict,dict) or not isinstance(kmer_comprev_query_dict,dict):
    raise ErroriPersonalizzati.NotADict()


sub = Subject('/Users/pietrodispaldro/Documents/GitHub/AlphaTeam/Esame/ref.fa')
diz_partenza_subject = sub.parse_file()
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
    '''
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
                                        '''
                        

    return seed_dict

a = find_seed(kmer_query_dict,kmer_subject_dict,kmer_comprev_subject_dict)
print(a)

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
#print(a)
#print(len(a))


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

#def extend_seed(schema,query.diz,):

#    return

transizione = {'A':'G','G':'A','C':'T','T':'C'}
trasversione = {'A':'C','A':'T','C':'A','C':'G','G':'C','G':'T','T':'A','T':'G'}


def extend_seed_right(schema,diz_partenza_query,diz_partenza_subject,k,x_max):
    score = 0
    match_consecutivi = 0

    for i in range(0,len(schema),3):
        prova = schema[i:i+3]
        query = prova[1][0]
        subject = prova[2][0]
        pos_query = prova[1][1]
        pos_sub = prova[2][1]
        hsp = prova[0]
        start_query = pos_query[0]
        start_subject = pos_sub[0]

        for j in diz_partenza_query.keys():
            if query == j:
                sequence_query = diz_partenza_query[j]

        for z in diz_partenza_subject.keys():
            if subject == z:
                sequence_sub = diz_partenza_subject[z]


        '''
        estensione_sinistra = ""
        mismatch_consecutivi = 0

        sequence_query_left = sequence_query[:start_query-k]
        sequence_subject_left = sequence_sub[:start_subject-k]

    for b in range(0,len(sequence_query_left)):
        if sequence_query_left[b] == sequence_subject_left[b]:
            estensione_sinistra = sequence_query_left[b] + estensione_sinistra
            mismatch_consecutivi = 0
            score +=1
        else:
            mismatch_consecutivi += 1
            chiave = sequence_query_left[b]
            if transizione[chiave] == sequence_subject_left[b]:
                score -= 1
            elif trasversione[chiave] == sequence_subject_left[b]:
                score -= 1
            if mismatch_consecutivi == x_max:
                print(f"Estensione a sinistra: Mi sono fermato in posizione : {b}")
                print(f"Sequenz: {sequence_query[0:b + 5]}")
                print(f"Subject: {sequence_sub[0:b + 5]}")
                break


        hsp = estensione_sinistra + hsp
        '''

        sequence_query = sequence_query[start_query + k:]
        sequence_sub = sequence_sub[start_subject + k:]


    for a in range(0, len(sequence_query)):
        if sequence_query[a] == sequence_sub[a]:
            match_consecutivi = 0
            score += 1
        else:
            match_consecutivi += 1
            chiave = sequence_query[a]
            if transizione[chiave] == sequence_sub[a]:
                score -= 1
            elif trasversione[chiave] == sequence_sub[a]:
                score -= 1
            if match_consecutivi == x_max:
                print(f"Estensione a destra: Mi sono fermato in posizione : {a}")
                print(f"Sequenz: {sequence_query[0:a-5]}")
                print(f"Subject: {sequence_sub[0:a-5]}")
                break

    hsp +=  sequence_query[0:a-(x_max-1)]
    score = score+k
    return sequence_query[0:a-(x_max-1)],score

def extend_seed_left(schema, diz_partenza_query, diz_partenza_subject, k, x_max):
        score = 0
        match_consecutivi = 0

        for i in range(0, len(schema), 3):
            prova = schema[i:i + 3]
            hsp = prova[0]
            query = prova[1][0]
            subject = prova[2][0]
            pos_query = prova[1][1]
            pos_sub = prova[2][1]
            start_query = pos_query[0]
            start_subject = pos_sub[0]


            for j in diz_partenza_query.keys():
                if query == j:
                    sequence_query = diz_partenza_query[j]

            for z in diz_partenza_subject.keys():
                if subject == z:
                    sequence_sub = diz_partenza_subject[z]


            sequence_query_left = sequence_query[:start_query]
            sequence_subject_left = sequence_sub[:start_subject]


            estensione_sinistra = ""
            mismatch_consecutivi = 0

            if start_query == 0:
                estensione_sinistra = ""
                print("Nessuna estensione possibile a sinistra")
            else:
                for queryBase, subjectBase in zip(
                        reversed(sequence_query_left),
                        reversed(sequence_subject_left)
                ):
                    if queryBase == subjectBase:
                        score += 1
                        mismatch_consecutivi = 0
                    else:
                        mismatch_consecutivi += 1
                        chiave = queryBase
                        if transizione[chiave] == subjectBase:
                            score -= 1
                        elif trasversione[chiave] == subjectBase:
                            score -= 1

                        if mismatch_consecutivi == x_max:
                            print(f"Estensione a sinistra:  Mi sono fermato dopo {mismatch_consecutivi} mismatch consecutivi.")
                            break


                    estensione_sinistra = queryBase + estensione_sinistra

                hsp = estensione_sinistra + hsp
                score += k

            return estensione_sinistra, score





print(extend_seed_right(schema,diz_partenza_query,diz_partenza_subject,22,6))
#print(extend_seed_left(schema,diz_partenza_query,diz_partenza_subject,22,6))



diz_partenza_query_def = diz_partenza_query
diz_partenza_subject_def = diz_partenza_subject

dizionario_seed = a
print(dizionario_seed)



schema = []
kmer = "TGAGGAATATTGGTCAATGGGC"
query_header = "b6635d67cb594473ddba9f8cfba5d13d"
subject_header = "MJ030-2-barcode67-umi101484bins-ubs-3"
pos_query = [0]
pos_subject = [330]

schema.append(kmer)
schema.append((query_header, pos_query))
schema.append((subject_header, pos_subject))

k = len(kmer)
x_max = 6


hsp_extended,score_extended = extend_seed_left(
    schema,
    diz_partenza_query_def,
    diz_partenza_subject_def,
    k,
    x_max
)

hsp_right,score_right = extend_seed_right(
    schema,
    diz_partenza_query,
    diz_partenza_subject,
    k,
    x_max
)
print("\nRisultato a sinistra")
print(kmer)
print(f"HSP Esteso senza seed: {hsp_extended}")
print(f"Score: {score_extended}")

print("\nRisultato a destra")
print(kmer)
print(f"HSP Esteso senza seed: {hsp_right}")
print(f"Score: {score_right}")

print("\n=== Risultato finale ===")
print(kmer)
print(hsp_extended+f'\033[1m{kmer}\033[0m'+hsp_right)
print(score_extended+score_right)









"""






#return

"""
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

