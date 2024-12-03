def hard_trimming(sequence:str, qualList: list[int], qs: int, len_s:int):
    """
    Ad un certo punto mettiamo la documentazione
    """
    # filter(lambda x,qs: x < qs, qualList)
    l = [x for x in qualList if x < qs]
    if l:
        i = qualList.index(l[0])
        qualList = qualList[:i]
        sequence = sequence[:i]
    if len(sequence) < len_s:
        sequence, qualList = None, None
    return sequence, qualList

def dynamic_trimming(seq: str, qscore: list[int], threshold_q: int = 20, min_length: float = 0.95) -> tuple[str, list[int]]:
    """Dynamic trimming using a sliding calculation."""
    len_seq = len(seq)
    if not len_seq:
        return "", []

    min_length_abs = round(len_seq*min_length)

    filter_list = [True] * len_seq
    i = 0

    while i < len_seq:
        if qscore[i] < threshold_q:
            filter_list[i] = False
            j = i + 1
            window_sum = qscore[i]
            window_size = 1

            while j < len_seq and (window_sum / window_size) < threshold_q:
                window_sum += qscore[j]
                window_size += 1
                filter_list[j] = False
                j += 1

            i = j
        else:
            i += 1

    if sum(filter_list) >= min_length_abs:
        return (''.join(s for s, f in zip(seq, filter_list) if f),
                [q for q, f in zip(qscore, filter_list) if f])
    return "", []