#!/usr/bin/env python


def complement(T):
    '''
    >>> complement('AAAACCCGGT')
    'ACCGGGTTTT'
    '''
    import string
    return T.translate(string.maketrans('ACGT', 'TGCA'))[::-1]

if __name__ == '__main__':
    import doctest
    doctest.testmod()
