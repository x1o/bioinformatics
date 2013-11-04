#!/usr/bin/env python

def find_approx_match_ix(P, S, d):
    '''
    >>> find_approx_match_ix('ATTCTGGA', 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 3)
    [6, 7, 26, 27]
    '''

    i = 0
    match_ix = []
    L = len(S)
    k = len(P)
    while i <= L - k:
        word = S[i:i+k]
        j = 0
        cur_mism_n = 0
        while (j <= k - 1) and (cur_mism_n <= d):
            if word[j] != P[j]:
                cur_mism_n += 1
            j += 1
        if cur_mism_n <= d:
            match_ix.append(i)
        i += 1

    return match_ix

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    print find_approx_match_ix('AAAAA', 'AACAAGCTGATAAACATTTAAAGAG', 1)
