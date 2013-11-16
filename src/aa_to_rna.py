#!/usr/bin/env python

def aa_to_rna(targ_aa):
    '''
    >>> aa_to_rna('MA')
    ['AUGGCA', 'AUGGCC', 'AUGGCG', 'AUGGCU']
    '''

    import itertools

    T = {}
    with open('RNA_codon_table_1.txt', 'r') as f:
        for line in f:
            try:
                codon, aa = line.strip().split()
                try:
                    T[aa].append(codon)
                except KeyError:
                    # First seen
                    T[aa] = [codon]
            except ValueError:
                # Hit a stop codon: ['UAA', 'UAG', 'UGA']
                pass

    cand_aax = [T[aa] for aa in targ_aa]
    return [''.join(pair) for pair in
            [tup for tup in itertools.product(*cand_aax)]]

if __name__ == '__main__':
    import doctest
    doctest.testmod()
