#!/usr/bin/env python

def num_of_peptides(target_mass):
    '''
    >>> num_of_peptides(1024)
    14712706211
    '''

    from gen_bv import gen_bin_vx
    import numpy as np

    mx = []
    with open('integer_mass_table.txt', 'r') as f:
        for line in f:
            mass = line.strip().split()[1]
            mx.append(int(mass))

    mx = [10, 2, 12]

    dot_p = lambda a, b: sum([p[0] * p[1] for p in zip(a, b)])
    sum_v = lambda a, b: [sum(p) for p in zip(a, b)]

    P = [[0 for _ in range(len(mx))]]  # permutations
    i = 0
    for i in range(1, len(mx)+1):
        P += gen_bin_vx([0 for _ in range(len(mx))], 0, i)

    # P = [np.array(v) for v in P]

    def gen_C(C, lvl, M=[]):
        print 'lvl %s' % lvl
        if lvl == 0:
            return M
        else:
            C_new = []
            for v in C[:]:
                for v_p in P:
                    s = sum_v(v, v_p)
                    # s = v + v_p
                    # print '%s + %s -> %s' % (v, v_p, s)
                    if s not in C_new:
                        if dot_p(v, mx) == target_mass:
                        # if np.dot(v, mx) == target_mass:
                            if v not in M:
                                M.append(v)
                        C_new.append(s)
            return gen_C(C_new, lvl-1, M)

    max_lvl = target_mass / min(mx)
    print 'max level: %s' % max_lvl

    C = P[:]
    from pprint import pprint
    pprint(sorted(gen_C(C, max_lvl)))
    # pprint(sorted(gen_C(C, 15)))

    return sorted(C)


if __name__ == '__main__':
    # import doctest
    # doctest.testmod()
    (num_of_peptides(50))
