# coding: utf-8

from base_map import BaseMap

def chain2proto(chain_code):
    """Convert chain code like 'ijK' to curve prototype."""
    
    dim = len(set(''.join(chain_code).lower()))

    assert dim <= 6
    letters = 'ijklmn'

    vect_dict = {}
    for k in range(dim):
        coord = [0]*dim
        coord[k] = 1
        vect_dict[letters[k]] = coord
        vect_dict[letters[k].upper()] = [-m for m in coord]
        
    def diag_coord(vector):
        arg = [vect_dict[k] for k in vector]
        coord = list(map(sum,zip(*arg)))
        return coord

    proto = [list(map(vect_dict.get,chain_code)) if len(chain_code) == 1 else diag_coord(m) for m in chain_code]
    
    proto = [[0] * dim] + proto
    for l in range(len(proto)-1):
        proto[l+1] = [c + d for c, d in zip(proto[l], proto[l+1])]

    return proto

def basis2base_map(basis):
    
    dim = len(basis)-1 if basis[-1] in ['0','1'] else len(basis)
    
    letters = 'ijklmn'
    assert dim <= 6

    l2i = {l: i for i, l in enumerate(letters)}
    perm = [None]*dim
    flip = [None]*dim
    time_rev = True if basis[-1] =='1' else False
    
    basis = basis[:-1] if basis[-1] in ['0','1'] else basis

    for k, l in enumerate(basis):
        lk = l.lower()
        perm[k] = l2i[lk]
        flip[k] = (l != lk)

    return BaseMap(perm, flip, time_rev)
