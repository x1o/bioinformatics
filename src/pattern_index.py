#!/usr/bin/env python

from complement import complement

def pat_ix(P, S):
    '''
    >>> pat_ix('ATAT', 'GATATATGCATATACTT')
    (1, 3, 9)
    >>> pat_ix('z', 'GATATATGCATATACTT')
    ()
    >>> pat_ix('AT', 'A')
    ()
    '''
    i = 0
    ix = []
    try:
        while True:
            ix.append(S.index(P, i))
            i = ix[-1] + 1
    except ValueError:
        return tuple(ix)

if __name__ == '__main__':
    # import doctest
    # doctest.testmod()
    # S = open('/home/xio/Desktop/Vibrio_cholerae.txt', 'r').read().strip()
    S = open('/home/xio/Desktop/Thermotoga-petrophila.txt', 'r').read().strip()
    P = 'CTTGATCAT'
    print ' '.join([str(el) for el in pat_ix(P, S)])
    print ' '.join([str(el) for el in pat_ix(complement(P), S)])
