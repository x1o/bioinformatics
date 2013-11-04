#!/usr/bin/env python

import itertools
import logging
dbg = logging.debug

def gen_mutations(S, d):
    '''
    >>> gen_mutations('AAC', 2)
    set(['ACC', 'ATG', 'AAG', 'AAA', 'ATC', 'AAC', 'ATA', 'AGG', 'ACT', 'AGC', 'ACA', 'AGA', 'CAT', 'AAT', 'ATT', 'CTC', 'CAC', 'ACG', 'CAA', 'AGT', 'CAG', 'CCC', 'TAT', 'CGC', 'GAT', 'TGC', 'TAG', 'TAA', 'GGC', 'TAC', 'TTC', 'GAC', 'TCC', 'GAA', 'GCC', 'GTC', 'GAG'])
    '''

    def has_inv(l):
        def is_inv(p):
            return p[0] >= p[1]

        if len(l) < 2:
            return False
        elif len(l) == 2:
            return is_inv(l)
        else:
            return is_inv(l[:2]) or has_inv(l[1:])

    # List of indices at which to mutate the DNA
    ix = [p for p in itertools.product(*[range(len(S)) for _ in range(d)]) if not has_inv(p)]
    # [A,T,C,G]^d (Cartesian product)
    nx = [p for p in itertools.product(*[['A', 'T', 'C', 'G'] for _ in range(d)])]
    # List of resulting strings with mutations
    sx = set([S])

    dbg('ix: %s' % ix)
    dbg('nx: %s' % nx)
    for idx in ix:
        for nucl_tup in nx:
            S_new = list(S)
            for i in range(d):
                S_new[idx[i]] = nucl_tup[i]
            sx.add(''.join(S_new))

    dbg('sx: %s' % sx)
    return sx


def count_k_mers(T, k, threshold=False, strict_eq=False, d=0):
    '''
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)
    ['CATG', 'GCAT']
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 5)
    []
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 2)
    ['CATG', 'ATGA', 'TGCA', 'GCAT']
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, 2, True)
    ['ATGA', 'TGCA']
    >>> count_k_mers('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4, False, False, 1)
    ['GATG', 'ATGC', 'ATGT']
    >>> count_k_mers('CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC', 10, False, False, 2)
    ['GCACACAGAC', 'GCGCACACAC']
    >>> count_k_mers('CTTGCCGGCGCCGATTATACGATCGCGGCCGCTTGCCTTCTTTATAATGCATCGGCGCCGCGATCTTGCTATATACGTACGCTTCGCTTGCATCTTGCGCGCATTACGTACTTATCGATTACTTATCTTCGATGCCGGCCGGCATATGCCGCTTTAGCATCGATCGATCGTACTTTACGCGTATAGCCGCTTCGCTTGCCGTACGCGATGCTAGCATATGCTAGCGCTAATTACTTAT', 9, False, False, 3)
    ['AGCGCCGCT', 'AGCGGCGCT']
    '''

    from operator import ge, eq

    i = 0
    L = len(T)
    D = {}
    dbg('>>> %s - %s - %s <<<' % (T, k, d))

    def add_kmer(D, kmer):
        if not D.has_key(kmer):
            dbg('> %s' % kmer)
            D[kmer] = 1
        else:
            dbg('+ %s' % kmer)
            D[kmer] += 1

    while i <= L-k:
        cur_p = T[i:i+k]
        dbg(' ~ %s' % cur_p)
        if d > 0:
            M = gen_mutations(cur_p, d)
            for kmer in M:
                add_kmer(D, kmer)
        else:
                add_kmer(D, cur_p)
        i += 1

    if not threshold:
        threshold = max(D.values())
    if strict_eq:
        op = eq
    else:
        op = ge

    dbg(sorted(D.items(), key=lambda el: el[1], reverse=True))
    return [item[0] for item in D.items() if op(item[1], threshold)]
    # return [item for item in D.items() if item[1] >= 3]

if __name__ == '__main__':
    logging.basicConfig(level='DEBUG')
    # logging.basicConfig(level='INFO')
    import doctest
    doctest.testmod()
    # print gen_mutations('ACTG', 1)
    # print gen_mutations('ACTG', 2)
    # S = 'AATGATGATGACGTCAAAAGGATCCGGATAAAACATGGTGATTGCCTCGCATAACGCGGTATGAAAATGGATTGAAGCCCGGGCCGTGGATTCTACTCAACTTTGTCGGCTTGAGAAAGACCTGGGATCCTGGGTATTAAAAAGAAGATCTATTTATTTAGAGATCTGTTCTATTGTGATCTCTTATTAGGATCGCACTGCCCTGTGGATAACAAGGATCCGGCTTTTAAGATCAACAACCTGGAAAGGATCATTAACTGTGAATGATCGGTGATCCTGGACCGTATAAGCTGGGATCAGAATGAGGGGTTATACACAACTCAAAAACTGAACAACAGTTGTTCTTTGGATAACTACCGGTTGATCCAAGCTTCCTGACAGAGTTATCCACAGTAGATCGCACGATCTGTATACTTATTTGAGTAAATTAACCCACGATCCCAGCCATTCTTCTGCCGGATCTTCCGGAATGTCGTGATCAAGAATGTTGATCTTCAGTG'
    # for i in range(1, 10):
        # print i
        # print ' ',  count_k_mers(S, i, 3)
    # print count_k_mers(S, 9, 2)
    # S = 'TGCTGCATTGATTGAATGGCATTGAATATGTGCATGAATGGCGGCTGCATTGGGCTGCATTGATGATTGATTGATTGGGCGGCAATAATTGCATTGAATAATATTGAATATGATGAATTGCTGCATGTGCGGCATGATTGATGGGCATTGTGCTGCATGGGCGGCATTGGGCTGCATGAATATTGAATATGGGCTGCGGCATGATGTGCTGCTGCATTGTGCGGCATTGAATGGCATGGGCGGCGGCGGCATGATGATGAATGGCTGCATTGAATATTGGGCAATATGAATGGCTGCATGATTGAATATGGGCTGCTGCAATAAT'
    # print count_k_mers(S, 9, False, False, 2)
    # S = 'AACATCCAGTCCTCTCCCTCTCGTCGTCCCGTACCTCCCGTCAATCCACCTCAGTCCCCATCCCCCCCCCCGTCCCTCGTGTGTGTTCCCCGTCTCCCCACCGTTCTCTCCCGTTCGTTCCCCAATCTCGTAGTTCCCGTTCGTCGTCGTAACCCCAGTTCGTCCATCCCCCCCCTCGTCGTGTATCACCTCCCCTCC'
    # print count_k_mers(S, 9, False, False, 2)

