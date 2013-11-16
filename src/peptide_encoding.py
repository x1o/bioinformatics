#!/usr/bin/env python

from aa_to_rna import aa_to_rna
from na_conv import rna_to_dna
from complement import rev_compl

def peptide_enc(dna, targ_aa):
    '''
    >>> peptide_enc('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA')
    ['ATGGCC', 'ATGGCC', 'GGCCAT']
    >>> with open('peptide_encoding_data.txt', 'r') as f:
    ...     _ = f.readline()
    ...     dna = f.readline().strip()
    ...     aa = f.readline().strip()
    >>> for rna_pat in sorted(peptide_enc(dna, aa)):
    ...     print rna_pat
    AAAGAAGTTTTCGAACCACATTATTAC
    AAAGAAGTTTTCGAGCCGCACTACTAC
    AAAGAGGTGTTTGAACCTCATTACTAT
    AAGGAAGTATTCGAACCACATTACTAT
    AAGGAAGTATTTGAGCCTCATTATTAC
    AAGGAAGTGTTTGAACCTCACTATTAT
    AAGGAGGTATTTGAACCCCACTATTAC
    ATAATAATGCGGCTCGAATACTTCCTT
    ATAATAATGTGGCTCGAACACTTCTTT
    ATAATAGTGAGGCTCAAAAACTTCCTT
    ATAGTAATGAGGTTCGAAAACCTCTTT
    ATAGTAATGGGGTTCGAAGACTTCCTT
    ATAGTAGTGAGGTTCGAAGACTTCCTT
    ATAGTAGTGTGGTTCAAATACCTCCTT
    GTAATAGTGCGGTTCAAAAACTTCCTT
    GTAGTAATGAGGTTCAAAAACCTCCTT
    GTAGTAATGGGGCTCAAACACCTCTTT
    GTAGTAATGGGGCTCGAAAACCTCCTT
    GTAGTAATGGGGTTCGAAGACTTCCTT
    GTAGTAGTGCGGCTCAAAAACTTCCTT
    '''

    def find_substrs(S, substr, sx):
        '''
        >>> find_substrs('abaaacaav', 'aa', [])
        ['aa', 'aa']
        '''

        ind = S.find(substr)
        if ind == -1:
            return sx
        else:
            sx.append(S[ind:ind+len(substr)])
            find_substrs(S[ind+len(substr):], substr, sx)
            return sx

    found_dnax = []

    for rna_pat in aa_to_rna(targ_aa):
        dna_pat = rna_to_dna(rna_pat)
        rev_compl_dna_pat = rev_compl(dna_pat)
        # print '> %s...' % dna_pat
        for found_dna_pat in find_substrs(dna, dna_pat, []):
            found_dnax.append(found_dna_pat)
        for found_rev_dna_pat in find_substrs(dna, rev_compl_dna_pat, []):
            found_dnax.append(found_rev_dna_pat)
        # XXX: naming could be better
        print len(found_dnax)

    return found_dnax

if __name__ == '__main__':
    # import doctest
    # doctest.testmod()
    # with open('/home/xio/Desktop/dataset_18_6(2).txt', 'r') as f:
        # dna = f.readline().strip()
        # aa = f.readline().strip()
    # for dna_pat in peptide_enc(dna, aa):
        # print dna_pat
    dna = ''
    with open('../dna/B_brevis.txt', 'r') as f:
        for line in f:
            dna += line.strip()
    print 'Done reading DNA with length %s' % len(dna)
    print len(peptide_enc(dna, 'VKLFPWFNQY'))
