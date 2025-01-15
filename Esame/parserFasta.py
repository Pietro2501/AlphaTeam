def parse_fasta(filepath):

    records = {}
    current_header = None
    current_sequence_lines = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            if not line:
                continue


            if line.startswith(">"):

                if current_header is not None:
                    records[current_header] = "".join(current_sequence_lines)
                    current_sequence_lines = []


                current_header = line[1:].strip()

            else:
                        
                current_sequence_lines.append(line.upper())
                        
                if current_header is not None:
                        records[current_header] = "".join(current_sequence_lines)

                    
        header=[]
        for current_header,current_sequence_lines in records.items():
            if "N" in current_sequence_lines:
                header.append(current_header)
        for current_header in header:
            del records[current_header]


        
    return records 


a=parse_fasta('C:\\Users\\rossa\\OneDrive\\Documenti\\GitHub\\AlphaTeam\\Esame\\ref.fa')
print(a)