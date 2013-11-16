#!/usr/bin/env python

def rna_to_aa(rna):
    '''
    >>> rna_to_aa('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
    'MAMAPRTEINSTRING'
    '''
    T = {}
    with open('RNA_codon_table_1.txt', 'r') as f:
        for line in f:
            try:
                codon, aa = line.strip().split()
                T[codon] = aa
            except ValueError:
                # Hit a stop codon: ['UAA', 'UAG', 'UGA']
                T[line.strip()] = ''

    def chunk(S, l):
        for i in xrange(0, len(S), l):
            yield S[i:(i+l)]

    return ''.join([T[codon] for codon in chunk(rna, 3)])

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    S = open('/home/xio/Desktop/dataset_18_3.txt', 'r').read().strip()
    print rna_to_aa(S)
