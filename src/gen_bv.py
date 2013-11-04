#!/usr/bin/env python

'''
>>> gen_bin_vx([0, 0, 0, 0, 0], 0, 1)
[[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]
>>> gen_bin_vx([0, 0, 0, 0, 0], 0, 3)
[[1, 1, 1, 0, 0], [1, 1, 0, 1, 0], [1, 1, 0, 0, 1], [1, 0, 1, 1, 0], [1, 0, 1, 0, 1], [1, 0, 0, 1, 1], [0, 1, 1, 1, 0], [0, 1, 1, 0, 1], [0, 1, 0, 1, 1], [0, 0, 1, 1, 1]]
'''

def gen_bin_vx(V, init_i, lvl):
    i = init_i
    vx = []

    while i <= len(V) - lvl:
        V_new = V[:]
        V_new[i] = 1
        if lvl == 1:
            vx.append(V_new)
        else:
            vx += gen_bin_vx(V_new, i+1, lvl-1)
        i += 1

    return vx

if __name__ == '__main__':
    import doctest
    doctest.testmod()
