#!/usr/bin/env python

def peptide_to_masses(peptide):
    '''
    >>> peptide_to_masses('GYM')
    [57, 163, 131]
    '''

    D = {}  # integer mass table
    with open('integer_mass_table.txt', 'r') as f:
        for line in f:
            aa, mass = line.strip().split()
            D[aa] = int(mass)

    return [D[aa] for aa in list(peptide)]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
