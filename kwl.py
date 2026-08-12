"""2-WL coherent closure, thin-Schur (TS) filtration, CFI construction. Exact integer arithmetic."""
import numpy as np
from itertools import combinations

def wl2(A):
    """Stable global 2-WL partition of V x V. Returns integer label matrix."""
    n = A.shape[0]
    # atomic: 0 diag, 1 arc, 2 other
    col = np.full((n, n), 2, dtype=np.int64)
    col[A > 0] = 1
    np.fill_diagonal(col, 0)
    while True:
        k = col.max() + 1
        # signature: for each (x,y): col(x,y) + multiset over z of (col(x,z),col(z,y))
        # encode pair colours as col(x,z)*k + col(z,y); count occurrences via matrix trick
        sigs = {}
        new = np.empty((n, n), dtype=np.int64)
        # build count tensors: for each pair (x,y), histogram over z of (col[x,z], col[z,y])
        # do it row by row to save memory
        for x in range(n):
            cxz = col[x]  # length n
            for y in range(n):
                czy = col[:, y]
                pairkey = cxz * k + czy
                pairkey_sorted = np.sort(pairkey)
                key = (col[x, y], pairkey_sorted.tobytes())
                if key not in sigs:
                    sigs[key] = len(sigs)
                new[x, y] = sigs[key]
        if len(sigs) == col.max() + 1 and np.array_equal(_canon(new), _canon(col)):
            return _canon(col)
        col = new

def _canon(col):
    """Canonical relabel by first occurrence."""
    n = col.shape[0]
    flat = col.ravel()
    _, idx = np.unique(flat, return_index=True)
    order = flat[np.sort(idx)]
    remap = {v: i for i, v in enumerate(order)}
    out = np.vectorize(remap.get)(col)
    return out.astype(np.int64)

def refine_by(col, mats):
    """Refine partition (label matrix col) by level sets of each matrix in mats and by transpose."""
    n = col.shape[0]
    keys = [col, col.T] + mats
    stacked = np.stack([m.astype(np.int64) for m in keys], axis=-1)  # n x n x k
    flat = stacked.reshape(n * n, -1)
    _, inv = np.unique(flat, axis=0, return_inverse=True)
    return _canon(inv.reshape(n, n))

def ts_partition(A, max_rounds=500):
    """Stable TS partition: start from {diag, arcs, rest}; alternate
    H (level sets + transpose) with M (multiply by current arc-supported blocks)."""
    n = A.shape[0]
    col = np.full((n, n), 2, dtype=np.int64)
    col[A > 0] = 1
    np.fill_diagonal(col, 0)
    col = _canon(col)
    arcmask = (A > 0)
    for _ in range(max_rounds):
        K = col.max() + 1
        blocks = [(col == c).astype(np.int64) for c in range(K)]
        arcblocks = [B for B in blocks if ((B > 0) <= arcmask).all()]
        prods = []
        for B in arcblocks:
            for C in blocks:
                prods.append(B @ C)
                prods.append(C @ B)
        new = refine_by(col, prods)
        if new.max() == col.max():
            return col
        col = new
    raise RuntimeError("TS did not stabilise")

def cfi(edges, nverts):
    """Even-subset CFI gadget graph. Vertices: ('E',e,b) and ('V',v,S frozenset).
    (v,S) ~ (e,b) iff v in e and b == [e in S]."""
    E = [tuple(sorted(e)) for e in edges]
    inc = {v: [e for e in E if v in e] for v in range(nverts)}
    verts = []
    for e in E:
        verts.append(('E', e, 0)); verts.append(('E', e, 1))
    for v in range(nverts):
        d = len(inc[v])
        for r in range(0, d + 1, 2):
            for S in combinations(inc[v], r):
                verts.append(('V', v, frozenset(S)))
    idx = {w: i for i, w in enumerate(verts)}
    n = len(verts)
    A = np.zeros((n, n), dtype=np.int64)
    for w in verts:
        if w[0] == 'V':
            v, S = w[1], w[2]
            for e in inc[v]:
                b = 1 if e in S else 0
                u = ('E', e, b)
                A[idx[w], idx[u]] = 1
                A[idx[u], idx[w]] = 1
    return A, verts

def analyze(A):
    """Return (K, dimTS, arc_rigid, diag_rigid): whether TS partition matches
    coherent partition on arcs and on the diagonal."""
    g = wl2(A)
    l = ts_partition(A)
    K = g.max() + 1
    dTS = l.max() + 1
    arcmask = (A > 0)
    # arc rigidity: on arc positions, does l-partition refine g-partition (finer or equal)?
    # l is coarser overall; agreement on arcs means: g-labels constant on l-classes restricted to arcs
    def agrees(mask):
        pairs = {}
        for lab_l, lab_g in zip(l[mask], g[mask]):
            if lab_l in pairs and pairs[lab_l] != lab_g:
                return False
            pairs[lab_l] = lab_g
        return True
    diagmask = np.eye(A.shape[0], dtype=bool)
    return K, dTS, agrees(arcmask), agrees(diagmask)
