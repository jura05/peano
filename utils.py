# coding: utf-8

from base_map import BaseMap

def chain2proto(chain_code, bases, start=None):
    """Convert chain code like 'ijK' to curve prototype."""
    
    if bases[0][-1] in ['0','1']:
        dim = len(bases[0])-1
    else:
        dim = len(bases[0])
    
    assert dim <= 6
    letters = 'ijklmn'
    l2v = {}
    for k in range(dim):
        l = letters[k]
        v = [0] * dim
        v[k] = 1
        l2v[l] = v
        l2v[l.upper()] = [-x for x in v]

    if start is None:
        start = (0,) * dim

    cube = start
    proto = [cube]
    for l in chain_code:
        diff = l2v[l]
        cube = [c + d for c, d in zip(cube, diff)]
        proto.append(cube)

    return proto

def basis2base_map(basis):
    dim = len(basis)
    letters = 'ijklmn'
    assert dim <= 6

    l2i = {l: i for i, l in enumerate(letters)}
    perm = [None]*dim
    flip = [None]*dim

    for k, l in enumerate(basis):
        lk = l.lower()
        perm[k] = l2i[lk]
        flip[k] = (l != lk)

    return BaseMap(perm, flip)
