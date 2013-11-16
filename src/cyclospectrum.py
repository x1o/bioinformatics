#!/usr/bin/env python

import logging
dbg = logging.debug

def cyclospectrum(peptide, p_type='cyclic'):
    '''
    >>> cyclospectrum('LEQN')
    [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
    >>> cyclospectrum('LEQN', p_type='linear')
    [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    '''

    A = [[], peptide]    # list of subpeptides as strings (e.g. ['GAS', 'LISP', ...])
    for ws in range(1, len(peptide)):   # window size
        i = 0
        while i <= len(peptide)-ws:
            dbg('%s:%s -> %s' % (i, i+ws, peptide[i:i+ws]))
            A.append(peptide[i:i+ws])
            i += 1
        if p_type == 'cyclic':
            for i in range(1, ws):
                dbg('%s + %s' % (peptide[-(ws-i):], peptide[:i]))
                A.append(peptide[-(ws-i):] + peptide[:i])
        dbg('')

    D = {}  # integer mass table
    with open('integer_mass_table.txt', 'r') as f:
        for line in f:
            aa, mass = line.strip().split()
            D[aa] = int(mass)

    pept_mass = lambda pept: sum([D[aa] for aa in list(pept)])

    return sorted([pept_mass(p) for p in A])


if __name__ == '__main__':
    # logging.basicConfig(level='DEBUG')
    logging.basicConfig(level='INFO')
    import doctest
    doctest.testmod()
    # for p in cyclospectrum('DTRTEYRYVHYRDAG'):
        # print p,
    print cyclospectrum('NAY')
