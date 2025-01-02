# Matrice per scoring transizioni e trasversioni
# trasversioni 3 volte più influenti delle transizioni
nucleotides_scoring_matrix = {
    ('A', 'A'): 1, ('A', 'C'): -3, ('A', 'G'): -1, ('A', 'T'): -3,
    ('C', 'A'): -3, ('C', 'C'): 1, ('C', 'G'): -3, ('C', 'T'): -1,
    ('G', 'A'): -1, ('G', 'C'): -3, ('G', 'G'): 1, ('G', 'T'): -3,
    ('T', 'A'): -3, ('T', 'C'): -1, ('T', 'G'): -3, ('T', 'T'): 1
}
