#!/usr/bin/env python


def rev_compl(T):
    '''
    >>> rev_compl('AAAACCCGGT')
    'ACCGGGTTTT'
    '''
    import string
    return T.translate(string.maketrans('ACGT', 'TGCA'))[::-1]

if __name__ == '__main__':
    import doctest
    doctest.testmod()
