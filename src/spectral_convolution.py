#!/usr/bin/env python

def spectral_convolution(spectrum):
    '''
    >>> spectral_convolution('0 137 186 323')
    [49, 137, 137, 186, 186, 323]
    '''
    spectrum = sorted([int(m) for m in spectrum.split(' ')])
    C = []
    i = 0
    while i < len(spectrum)-1:
        j = i+1
        while j <= len(spectrum)-1:
            C.append(spectrum[j] - spectrum[i])
            j += 1
        i += 1
    return sorted(filter(lambda el: el != 0, C))

def spectral_convolution_2(spectrum):
    spectrum = sorted([int(m) for m in spectrum.split(' ')])
    from itertools import product
    return filter(lambda m: m > 0, [p[0] - p[1] for p in product(spectrum, spectrum)])

def freq(convolution):
    M = set(convolution)
    D = {}
    for m in M:
        if m >= 50 and m <= 200:
            D[m] = convolution.count(m)
    return D

if __name__ == '__main__':
    # import doctest
    # doctest.testmod()
    # print ' '.join([str(m) for m in spectral_convolution('935 201 490 289 250 87 761 962 636 300 950 771 874 597 387 1049 413 1021 436 71 460 333 490 347 832 626 436 700 200 950 886 1136 1065 262 676 847 186 749 113 613 97 836 510 114 313 533 500 539 946 115 190 147 646 261 700 803 823 740 365 760 603 210 376 1022 174 674 304 646 1023 299 926 462 499 375 1039 637 1033 523 837 561 103 575 103 396 1033 789 186 875 989 0 723 936')])
    # print sorted(spectral_convolution_2('935 201 490 289 250 87 761 962 636 300 950 771 874 597 387 1049 413 1021 436 71 460 333 490 347 832 626 436 700 200 950 886 1136 1065 262 676 847 186 749 113 613 97 836 510 114 313 533 500 539 946 115 190 147 646 261 700 803 823 740 365 760 603 210 376 1022 174 674 304 646 1023 299 926 462 499 375 1039 637 1033 523 837 561 103 575 103 396 1033 789 186 875 989 0 723 936'))
    # print sorted(freq(spectral_convolution('0 113 114 128 129 227 242 242 257 355 356 370 371 484')).items(), key=lambda item: item[1], reverse=True)
    print ' '.join([str(m) for m in spectral_convolution('199 500 97 429 113 244 454 131 228 323 726 482 210 301 695 341 369 595 639 113 0 394 312 298 752 184 71 256 185 624 579 710 567 341 710 170 57 638 613 369 128 425 525 710 653 597 226 398 692 113 526 454 522 297 823 766 511 482')])
