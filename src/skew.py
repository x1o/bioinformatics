#!/usr/bin/env python

# Simple and beautiful O(n)

def skew(S):
    '''
    >>> skew('CATGGGCATCGGCCATACGCC')
    [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]
    '''

    nx = [0]
    for nucl in S:
        if nucl == 'C':
            nx.append(nx[-1] - 1)
        elif nucl == 'G':
            nx.append(nx[-1] + 1)
        else:
            nx.append(nx[-1])
        # sys.stdout.write('\r%.1f%%' % (float(i) / len(S) * 100))

    return nx

def plot_skew(nx):
    from matplotlib import pylab
    pylab.plot(nx)
    pylab.show()

def find_min_i(nx):
    '''
    >>> find_min_i([0, 1, 1, 0, 0, 1, 0, -1, 0, 0, -1, 1, 2, 3, 2, 1, 1, 1, 0, 0])
    [7, 10]
    '''

    min_ix = [0]
    i = 1
    while i <= len(nx) - 1:
        if nx[i] < nx[min_ix[-1]]:
            min_ix = [i]
        elif nx[i] == nx[min_ix[-1]]:
            min_ix.append(i)
        i += 1

    return min_ix


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # S = open('../dna/E-coli.txt', 'r').read().strip()
    # print ' '.join([str(n) for n in skew('GAGCCACCGCGATA')])
    # S = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
    # S = open('/home/xio/Desktop/dataset_7_6.txt', 'r').read().strip()
    # S = open('../dna/E-coli.txt', 'r').read().strip()
    # print ' '.join([str(n) for n in find_min_i(skew(S))])
    # print skew(S)
    # plot_skew(skew(S))
