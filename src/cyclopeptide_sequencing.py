#!/usr/bin/env python

from cyclospectrum import cyclospectrum
from peptide_to_masses import peptide_to_masses
import logging
dbg = logging.debug
info = logging.info

A = []  # AA's
with open('integer_mass_table.txt', 'r') as f:
    for line in f:
        aa = line.strip().split()[0]
        A.append(aa)

def expand(L):
    '''
    >>> expand({'': True}).keys()
    ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
    '''

    # k = len(L)
    # while k > 0:
        # peptide = L[0]
        # for aa in A:
            # L.append(peptide+aa)
        # L.remove(peptide)
        # k -= 1

    # XXX: Damn lists are slow
    # L_n = []
    # for peptide in L:
        # for aa in A:
            # L_n.append(peptide+aa)

    L_n = {}
    for peptide in L.iterkeys():
        for aa in A:
            L_n[peptide+aa] = -1

    return L_n

def is_consistent(peptide, spectrum):
    '''
    (HAY vs NAY)
    >>> is_consistent('HA', [0, 71, 137, 163, 208, 234, 300, 371])
    True
    >>> is_consistent('HA', [0, 71, 114, 163, 185, 234, 277, 348])
    False
    '''

    spectrum = spectrum[:]
    for mass in cyclospectrum(peptide, p_type='linear'):
        if mass not in spectrum:
            return False
        else:
            spectrum.remove(mass)   # XXX: should slow down the algorithm
    return True

def cyclopeptide_sequencing(spectrum):
    '''
    >>> cyclopeptide_sequencing('0 113 128 186 241 299 314 427')
    [[113, 128, 186], [113, 186, 128], [128, 113, 186], [128, 186, 113], [186, 113, 128], [186, 128, 113]]
    '''

    spectrum = [int(mass) for mass in spectrum.split()]
    # L = ['']

    L = {}
    L[''] = -1
    M = []  # matches

    while L:
        L = expand(L)
        dbg([peptide_to_masses(peptide) for peptide in L.iterkeys()])
        # print [peptide_to_masses(peptide) for peptide in L]
        dbg('* %s' % spectrum)
        k = len(L)
        dbg('List of length %s' % k)

        # for peptide in L:
        for peptide in L.keys():
            if cyclospectrum(peptide) == spectrum:
                M.append(peptide)
                info('*** Found %s' % peptide)
                # L.remove(peptide)
                del L[peptide]
            elif not is_consistent(peptide, spectrum):
                # L.remove(peptide)
                del L[peptide]
        dbg('Removed %s %% peptides' % ((k - len(L)) / float(k) * 100))
        dbg(sorted([peptide_to_masses(peptide) for peptide in L.iterkeys()]))
        # print sorted([peptide_to_masses(peptide) for peptide in L])
        dbg('')

    # filter out duplicated matches
    M_masses = []
    for peptide in M:
        masses = peptide_to_masses(peptide)
        if masses not in M_masses:
            M_masses.append(masses)

    return sorted(M_masses)

def leaderboard_cyclopeptide_sequencing(spectrum, N):
    '''
    >>> leaderboard_cyclopeptide_sequencing('0 71 113 129 147 200 218 260 313 331 347 389 460', 10)
    [[113, 147, 71, 129]]
    '''

    def score(peptide, spectrum):
        '''
        >>> spectrum = '0 97 99 114 128 147 147 163 186 227 241 242 244 260 261 262 283 291 333 340 357 385 389 390 390 405 430 430 447 485 487 503 504 518 543 544 552 575 577 584 632 650 651 671 672 690 691 738 745 747 770 778 779 804 818 819 820 835 837 875 892 917 932 932 933 934 965 982 989 1030 1039 1060 1061 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1225 1322'
        >>> score('VKLFPWFNQY', [int(m) for m in spectrum.split()])
        86
        '''

        # a little cryptic; should be fast though
        i = j = 0
        res = 0   # score
        es = cyclospectrum(peptide) # experimental spectrum
        ts = spectrum    # theoretical spectrum
        while (i <= len(ts)-1) and (j <= len(es)-1):
            if ts[i] == es[j]:
                res += 1
                i += 1
                j += 1
            elif ts[i] > es[j]:
                j += 1
            elif ts[i] < es[j]:
                i += 1
        return res

    def cut(L, N):
        '''
        >>> cut({'a': 1, 'b': 2, 'c': 2, 'd': 3}, 2)
        {'c': 2, 'b': 2, 'd': 3}
        '''
        L = sorted(L.items(), key = lambda p: p[1], reverse=True)
        i = N-1
        try:
            while (L[i][1] == L[i+1][1]) and (L[i][1] != -1):
                i += 1
        except IndexError as e:
            # print 'IndexError: %s' % e
            pass
        # print 'Reducing to the first %s elements out of %s' % (i+1, len(L))
        # print L[0:i+1]
        return dict(L[0:i+1])


    mass = lambda peptide: sum(peptide_to_masses(peptide))
    parent_mass = lambda spectrum: sorted(spectrum)[-1]

    spectrum = [int(m) for m in spectrum.split()]

    L = {}  # leaderboard
    L[''] = -1    # add 0-peptide
    leader = ''     # current leader peptide

    while L:
        L = expand(L)
        print 'List of length %s' % len(L)

        for peptide in L.keys():
            p_score = score(peptide, spectrum)
            L[peptide] = p_score
            # if p_score >= 83:
                # print '%s (%s)' % (peptide, p_score)
            if mass(peptide) == parent_mass(spectrum):
                if p_score > score(leader, spectrum):
                    print 'Leader -> %s (%s)' % (peptide, score(peptide, spectrum))
                    leader = peptide
            elif mass(peptide) > parent_mass(spectrum):
                del L[peptide]

        # print 'Cutting...'
        # print sorted(L.items(), key=lambda p: p[1], reverse=True)
        # print '|'
        # print 'V'
        L = cut(L, N)
        print 'Current leader: ', leader
        if L:
            W = sorted(L.items(), key=lambda p: p[1], reverse=True)
            print '20 highest scorers: ', W[0:20]
        print

    print 'L is empty'
    print spectrum
    # print cyclospectrum(w[0])
    print leader
    return peptide_to_masses(leader)


if __name__ == '__main__':
    # logging.basicConfig(level='DEBUG', format='')
    logging.basicConfig(level='INFO')
    # import doctest
    # doctest.testmod()
    # print cyclopeptide_sequencing('0 97 97 99 101 103 196 198 198 200 202 295 297 299 299 301 394 396 398 400 400 497')
    # print '-'.join([str(m) for m in leaderboard_cyclopeptide_sequencing('0 71 113 129 147 200 218 260 313 331 347 389 460', 10)])
    # print '-'.join([str(m) for m in leaderboard_cyclopeptide_sequencing('0 87 87 87 97 99 99 99 113 113 113 115 128 128 128 129 131 137 137 156 156 163 202 202 210 215 216 218 225 227 236 236 241 241 242 243 255 255 262 265 268 269 276 289 328 329 330 333 338 340 342 355 362 364 366 367 368 370 372 375 390 392 392 404 417 420 443 454 457 461 465 469 470 471 475 479 485 491 494 503 503 503 505 520 523 530 548 556 557 569 574 584 590 593 602 606 607 608 610 613 616 617 631 631 632 636 656 661 685 693 704 705 712 712 718 719 721 723 725 730 730 730 730 733 736 744 745 792 798 811 812 817 820 824 829 831 832 833 833 838 849 849 861 867 873 875 886 897 920 923 925 935 940 946 946 946 948 948 960 960 961 966 973 974 977 985 988 998 1022 1033 1047 1051 1053 1053 1059 1060 1060 1063 1072 1074 1075 1076 1085 1087 1097 1102 1111 1116 1122 1150 1150 1159 1162 1162 1164 1166 1173 1184 1187 1188 1200 1201 1203 1209 1215 1215 1216 1224 1239 1249 1253 1263 1278 1286 1287 1287 1293 1299 1301 1302 1314 1315 1318 1329 1336 1338 1340 1340 1343 1352 1352 1380 1386 1391 1400 1405 1415 1417 1426 1427 1428 1430 1439 1442 1442 1443 1449 1449 1451 1455 1469 1480 1504 1514 1517 1525 1528 1529 1536 1541 1542 1542 1554 1554 1556 1556 1556 1562 1567 1577 1579 1582 1605 1616 1627 1629 1635 1641 1653 1653 1664 1669 1669 1670 1671 1673 1678 1682 1685 1690 1691 1704 1710 1718 1757 1758 1766 1769 1772 1772 1772 1772 1777 1779 1781 1783 1784 1790 1790 1797 1798 1809 1817 1841 1846 1866 1870 1871 1871 1885 1886 1889 1892 1894 1895 1896 1900 1909 1912 1918 1928 1933 1945 1946 1954 1972 1979 1982 1997 1999 1999 1999 2008 2011 2017 2023 2027 2031 2032 2033 2037 2041 2045 2048 2059 2082 2085 2098 2110 2110 2112 2127 2130 2132 2134 2135 2136 2138 2140 2147 2160 2162 2164 2169 2172 2173 2174 2213 2226 2233 2234 2237 2240 2247 2247 2259 2260 2261 2261 2266 2266 2275 2277 2284 2286 2287 2292 2300 2300 2339 2346 2346 2365 2365 2371 2373 2374 2374 2374 2387 2389 2389 2389 2403 2403 2403 2405 2415 2415 2415 2502', 381)])
    print leaderboard_cyclopeptide_sequencing('0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322', 1000)
