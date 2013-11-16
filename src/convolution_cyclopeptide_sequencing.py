#!/usr/bin/env python

from spectral_convolution import spectral_convolution, freq
from cyclopeptide_sequencing import leaderboard_cyclopeptide_sequencing

def convolution_cyclopeptide_sequencing(M, N, spectrum):
    '''
    >>> convolution_cyclopeptide_sequencing(20, 60, '57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493')
    [99, 71, 137, 57, 72, 57]
    '''
    # alphabet
    A = sorted(freq(spectral_convolution(spectrum)).items(), key=lambda it: it[1], reverse=True)
    i = M-1
    try:
        while (A[i][1] == A[i+1][1]):
            i += 1
    except IndexError:
        pass
    A = [p[0] for p in A[0:i+1]]

    # XXX
    # leaderboard_cyclopeptide_sequencing(...)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
