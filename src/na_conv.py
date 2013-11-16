#!/usr/bin/env python

import string

def rna_to_dna(rna):
    '''
    >>> rna_to_dna('GCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
    'GCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    '''
    return rna.translate(string.maketrans('U', 'T'))

def dna_to_rna(dna):
    '''
    >>> dna_to_rna('GCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA')
    'GCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    '''
    return dna.translate(string.maketrans('T', 'U'))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
