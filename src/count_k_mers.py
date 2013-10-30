#!/usr/bin/env python

def count_k_mers(T, k, threshold=False, exact=False):
    '''
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)
    ('CATG', 'GCAT')
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 5)
    ()
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 2)
    ('CATG', 'ATGA', 'TGCA', 'GCAT')
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 2, True)
    ('ATGA', 'TGCA')
    '''

    from operator import ge, eq

    i = 0
    L = len(T)
    D = {}
    while i < L-k:
        cur_word = T[i:i+k]
        if not D.has_key(cur_word):
            D[cur_word] = 1
        else:
            D[cur_word] += 1
        i += 1

    # while i < L-k:
        # cur_word = T[i:i+k]
        # # print '%s:' % cur_word
        # if not D.has_key(cur_word):
            # D[cur_word] = 1
            # j = i + 1
            # while j <= L-k and not D[cur_word] == max_t:
                # # print '  %s' % T[j:j+k]
                # if cur_word == T[j:j+k]:
                    # # print '  eq!'
                    # D[cur_word] += 1
                    # # print D
                # j += 1
        # i += 1

    if not threshold:
        threshold = max(D.values())
    if exact:
        op = eq
    else:
        op = ge
    # print sorted(D.items(), key=lambda el: el[1], reverse=True)
    return tuple([item[0] for item in D.items() if op(item[1], threshold)])
    # return tuple([item for item in D.items() if item[1] >= 3])

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # S = 'AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTTGTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTTAGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCGTTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAACCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA'
    # for i in range(1, 10):
        # print i
        # print ' ',  count_k_mers(S, i)
    # print count_k_mers(S, 9, 3)
