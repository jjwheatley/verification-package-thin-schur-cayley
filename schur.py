"""Group-level computation of the coherent closure and thin-Schur closure of a
connected Cayley graph over a finite abelian group.

Justification: both the coherent partition (lem:trans, Part III) and the stable
local partition (prop:trans, Part I) of a Cayley graph are translation
invariant, so each is determined by a partition of the group into "basic sets"
(difference classes).  The coherent closure's span is the smallest Schur ring
containing 1_S, computed by partition refinement under convolution of all class
pairs; TS is computed by the same refinement restricted to convolution by the
partition's own classes contained in S (the visible arc classes), plus
inversion.  Both agree with the O(n^2)-pair matrix computations (validated in
validate() below against kwl.wl2 / kwl.ts_partition).
"""
import numpy as np
import itertools

def _group(dims):
    els = list(itertools.product(*[range(d) for d in dims]))
    idx = {e: i for i, e in enumerate(els)}
    return els, idx

def _conv_table(dims):
    """index table: add[i,j] = index of els[i]+els[j]."""
    els, idx = _group(dims)
    n = len(els)
    add = np.zeros((n, n), dtype=np.int32)
    for i, x in enumerate(els):
        for j, y in enumerate(els):
            add[i, j] = idx[tuple((a + b) % d for a, b, d in zip(x, y, dims))]
    neg = np.array([idx[tuple((-a) % d for a, d in zip(x, dims))] for x in els],
                   dtype=np.int32)
    return els, idx, add, neg

def _canon(v):
    _, lab = np.unique(v, return_inverse=True)
    return lab

def _refine(lab, *vals):
    key = lab.astype(np.int64)
    for v in vals:
        key = key * (int(v.max()) + 1) + v
        _, key = np.unique(key, return_inverse=True)
    return key

def _conv_class(members_i, indicator, add):
    """convolution 1_{class} * vec at group level: out[k] = sum over i in class
    of indicator[k - i]  (uses precomputed addition table)."""
    out = np.zeros(len(indicator), dtype=np.int64)
    for i in members_i:
        out[add[i]] += indicator          # add[i] = row of sums i + (.)
    return out

def closures(dims, S, need_ts=True, add_tables=None):
    """Return (coh_labels, ts_labels) as arrays of class labels over the group
    elements (basic-set partitions).  S: iterable of group-element tuples,
    assumed inverse-closed; connectivity is NOT checked here."""
    if add_tables is None:
        els, idx, add, neg = _conv_table(dims)
    else:
        els, idx, add, neg = add_tables
    n = len(els)
    Sset = set(map(tuple, S))
    Sind = np.array([1 if e in Sset else 0 for e in els], dtype=np.int64)
    base = np.full(n, 2, dtype=np.int64)
    base[Sind > 0] = 1
    base[idx[tuple(0 for _ in dims)]] = 0

    def stabilise(lab, arcs_only):
        lab = _canon(lab)
        while True:
            K = lab.max() + 1
            classes = [np.where(lab == c)[0] for c in range(K)]
            if arcs_only:
                gen = [m for m in classes if Sind[m].all()]
            else:
                gen = classes
            C = np.zeros((K, n), dtype=np.int64)
            C[lab, np.arange(n)] = 1
            feats = [lab[None, :], lab[neg][None, :]]
            for m in gen:
                M = np.zeros((K, n), dtype=np.int64)
                for i in m:
                    M[:, add[i]] += C
                feats.append(M)
            F = np.vstack(feats)
            _, new = np.unique(F.T, axis=0, return_inverse=True)
            if new.max() == lab.max():
                return new
            lab = new

    coh = stabilise(base, arcs_only=False)
    ts = stabilise(base, arcs_only=True) if need_ts else None
    return coh, ts

def validate():
    """Cross-check against the matrix-level implementations on known cases."""
    from kwl import wl2, ts_partition
    def cay(dims, S):
        els, idx = _group(dims)
        n = len(els)
        A = np.zeros((n, n), dtype=np.int64)
        for i, x in enumerate(els):
            for s in S:
                A[i, idx[tuple((a + b) % d for a, b, d in zip(x, s, dims))]] = 1
        return A
    def pm(n, gens):   # circulant +-closure
        return [ (g % n,) for g in gens ] + [ ((-g) % n,) for g in gens ]
    cases = [
        ((24,), pm(24, [1, 5])),                     # C24(1,5): 9 / 9
        ((25,), pm(25, [1, 2])),                     # C25(1,2): 13 / 13
        ((24,), [(s,) for s in [2,3,10,14,21,22]]),  # wedge witness: 8 / 8
        ((13,), pm(13, [1, 5])),                     # prime: 4 / 4
        ((4, 8), [(1,0),(1,5),(1,6),(3,0),(3,2),(3,3)]),  # counterexample: 13 / 12
        ((2, 2, 6), [(1,0,0),(0,1,0),(0,0,1),(0,0,5),(1,1,3)]),
    ]
    for dims, S in cases:
        A = cay(dims, S)
        g = wl2(A); l = ts_partition(A)
        coh, ts = closures(dims, S)
        Km, Tm = g.max()+1, l.max()+1
        Kg, Tg = coh.max()+1, ts.max()+1
        # compare partitions of difference classes, not just counts
        row = g[0]; rowl = l[0]
        ok = (Km == Kg and Tm == Tg
              and len(set(zip(row.tolist(), coh.tolist()))) == Kg
              and len(set(zip(rowl.tolist(), ts.tolist()))) == Tg)
        print(f'{dims} |S|={len(S)}: matrix K={Km},T={Tm}  group K={Kg},T={Tg}  '
              f'{"OK" if ok else "MISMATCH"}')
        assert ok
    print('validation passed')

if __name__ == '__main__':
    validate()
