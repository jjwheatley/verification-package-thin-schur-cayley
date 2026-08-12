"""Verify every computational claim of "Thin-Schur closures of Cayley graphs:
positive circulant regimes, a minimal abelian counterexample, and its census".

Each claim prints a PASS or FAIL line; the script exits nonzero if any fails.
Runtime: ~4-5 minutes (add --full for the complete 995-graph atlas sweep).

    python3 verify.py [--full]
"""
import numpy as np, itertools, sys
sys.path.insert(0, '.')
from kwl import wl2, ts_partition, cfi
from certificates import classes, merged_pair, _rref_add, P1, P2
from schur import closures, _conv_table

ok = lambda b, msg: print(('PASS  ' if b else '*** FAIL  ') + msg) or b
allok = True

def circ(n, gens):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for g in gens:
            A[i, (i + g) % n] = 1; A[i, (i - g) % n] = 1
    return A

def circ_set(n, S):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for s in S: A[i, (i + s) % n] = 1
    return A

# ---- thm:prime spot checks: walk separation on prime circulants ----
for n, gens, Kexp in [(13, [1, 5], 4), (17, [1, 2, 4, 8], 3), (19, [1, 7, 8], 4)]:
    A = circ(n, gens); g = wl2(A); K = g.max() + 1
    P = np.eye(n, dtype=object); prof = [[] for _ in range(n)]
    for k in range(n):
        for y in range(n): prof[y].append(int(P[0, y]))
        P = P @ A
    reps = {}
    for y in range(n): reps.setdefault(int(g[0, y]), y)
    sep = len({tuple(prof[y]) for y in reps.values()}) == K
    allok &= ok(K == Kexp and sep, f'prime-order walk separation: C{n}{tuple(gens)}: K={K}, walk-separated')

# ---- C25(1,2): degenerate spectrum, walk separation holds ----
A = circ(25, [1, 2]); g = wl2(A); K = g.max() + 1
evs = {round(sum(2*np.cos(2*np.pi*a*s/25) for s in [1, 2]), 9) for a in range(25)}
P = np.eye(25, dtype=object); prof = [[] for _ in range(25)]
for k in range(25):
    for y in range(25): prof[y].append(int(P[0, y]))
    P = P @ A
reps = {}
for y in range(25): reps.setdefault(int(g[0, y]), y)
sep = len({tuple(prof[y]) for y in reps.values()}) == K
allok &= ok(K == 13 and len(evs) == 12 and sep,
            f'degenerate spectrum, C25(1,2): K={K}, {len(evs)} eigenvalues (degenerate), walk-separated')

# ---- prop:strict / prop:wa: the three wedge witnesses ----
for n, S, v1, v2, wts in [(24, [2,3,10,14,21,22], 33, 36, {15,18}),
                          (40, [1,5,9,11,19,21,29,31,35,39], 240, 222, {51,60}),
                          (48, [1,15,23,25,33,47], 36, 33, {15,18})]:
    A = circ_set(n, S); g = wl2(A); l = ts_partition(A); K = g.max()+1
    P = np.eye(n, dtype=object); prof = [[] for _ in range(n)]
    for k in range(n):
        for y in range(n): prof[y].append(int(P[0, y]))
        P = P @ A
    reps = {}
    for y in range(n): reps.setdefault(int(g[0, y]), y)
    byprof = {}
    for c, y in reps.items(): byprof.setdefault(tuple(prof[y]), []).append(y)
    ties = [v for v in byprof.values() if len(v) > 1]
    A3 = A@A@A; H = (A*A3)@A
    tie_ok = len(ties) == 1 and {int(H[0, ties[0][0]]), int(H[0, ties[0][1]])} == {v1, v2}
    arcs = {c: y for c, y in reps.items() if A[0, y]}
    w = {int((A*A3)[0, y]) for y in arcs.values()}
    allok &= ok(K == l.max()+1 and tie_ok and len(arcs) == 2 and w == wts,
                f'wedge witness n={n}: TS=Coh, one walk-tie, (AoA3)A={{{v1},{v2}}}, weights {sorted(wts)}')

# ---- lem:bmdual refutation on C24(1,5), exact cyclotomic arithmetic ----
n = 24
def red(k):
    poly = {k % 24: 1}; ch = True
    while ch:
        ch = False
        for d in sorted(poly, reverse=True):
            if d >= 8 and poly[d]:
                c = poly.pop(d); poly[d-4] = poly.get(d-4, 0)+c; poly[d-8] = poly.get(d-8, 0)-c; ch = True
    v = [0]*8
    for d, c in poly.items():
        if c: v[d] += c
    return tuple(v)
Z = [red(k) for k in range(24)]
def conv(u, v):
    w = np.zeros(24, dtype=object)
    for i in range(24):
        if u[i]: w = w + u[i]*np.roll(v, i)
    return w
def algdim(gens, p, unital=True):
    basis = []; vecs = []
    def add(v):
        vv = np.array([int(x) % p for x in v], dtype=np.int64)
        if _rref_add(basis, vv, p): vecs.append(v); return True
        return False
    if unital:
        e = np.zeros(24, dtype=object); e[0] = 1; add(e)
    for G in gens: add(G)
    i = 0
    while i < len(vecs):
        x = vecs[i]; i += 1
        for G in gens: add(conv(x, G))
    return len(vecs), vecs
arc = np.array([1 if y in (1,5,19,23) else 0 for y in range(24)], dtype=object)
dT, vecs = algdim([arc], P1)
A = circ(24, [1, 5]); g = wl2(A); K = g.max()+1
reps = {}
for y in range(24): reps.setdefault(int(g[0, y]), y)
profs = {c: tuple(int(v[y]) for v in vecs) for c, y in reps.items()}
sep = len(set(profs.values()))
def transform(v, b):
    s = np.zeros(8, dtype=object)
    for x in range(24):
        if v[x]: s = s + int(v[x])*np.array(Z[(b*x) % 24], dtype=object)
    return tuple(int(c) for c in s)
Dmap = {}
for b in range(24): Dmap.setdefault(tuple(transform(v, b) for v in vecs), []).append(b)
Delta = list(Dmap.values())
gens2 = [np.array([1 if b in blk else 0 for b in range(24)], dtype=object) for blk in Delta]
d1, _ = algdim(gens2, P1, unital=False); d2, _ = algdim(gens2, P2, unital=False)
allok &= ok(K == 9 and dT == 7 and sep == 9 and len(Delta) == 7 and d1 == d2 == 9,
            f'Bannai-Muzychuk cardinality gap on C24(1,5): dimT={dT}, |Lambda|={sep}, |Delta|={len(Delta)}, dual algebra dim={d1} (not a fusion)')

# ---- thm:counter: the Z4xZ8 counterexample, matrix level + certificate ----
dims = (4, 8); n = 32
els = list(itertools.product(range(4), range(8))); idx = {e: i for i, e in enumerate(els)}
def cay(S):
    A = np.zeros((n, n), dtype=np.int64)
    for i, x in enumerate(els):
        for s in S: A[i, idx[((x[0]+s[0]) % 4, (x[1]+s[1]) % 8)]] = 1
    return A
def fused_cert(A, g, c1, c2):
    f = g.copy(); f[f == c2] = c1
    _, f = np.unique(f, return_inverse=True); f = f.reshape(n, n)
    blocks = [(f == c).astype(np.int64) for c in range(f.max()+1)]
    def in_span(M):
        return all((M[B > 0] == M[B > 0][0]).all() for B in blocks if B.sum())
    arcs = [B for B in blocks if (B <= (A > 0)).all() and B.sum()]
    mult = all(in_span(R@B) and in_span(B@R) for R in arcs for B in blocks)
    tr = all(in_span(B.T) for B in blocks)
    iaj = in_span(np.eye(n, dtype=np.int64)) and in_span(A) and in_span(np.ones((n, n), dtype=np.int64))
    viol = sum(1 for B1 in blocks for B2 in blocks if not in_span(B1@B2))
    return mult and tr and iaj and viol > 0
S0 = [(1,0),(1,5),(1,6),(3,0),(3,2),(3,3)]
A = cay(S0); g = wl2(A); l = ts_partition(A)
c1 = int(g[0, idx[(0,2)]]); c2 = int(g[0, idx[(2,0)]])
merged = l[0, idx[(0,2)]] == l[0, idx[(2,0)]] and c1 != c2
allok &= ok(g.max()+1 == 13 and l.max()+1 == 12 and merged and fused_cert(A, g, c1, c2),
            'counterexample Z4xZ8 |S|=6: K=13, dimTS=12, merged pair certified (internal stall)')
# T fails to separate (group algebra, exact)
tabs = _conv_table(dims)
def gconv(u, v):
    w = np.zeros(n, dtype=object)
    for i in range(n):
        if u[i]:
            sh = np.zeros(n, dtype=object)
            for j in range(n):
                if v[j]: sh[tabs[2][i, j]] += v[j]
            w = w + u[i]*sh
    return w
R1 = np.array([1 if els[i] in [(1,0),(1,6),(3,0),(3,2)] else 0 for i in range(n)], dtype=object)
R2 = np.array([1 if els[i] in [(1,5),(3,3)] else 0 for i in range(n)], dtype=object)
basis = []; vecs = []
def add(v):
    vv = np.array([int(x) % P1 for x in v], dtype=np.int64)
    if _rref_add(basis, vv, P1): vecs.append(v); return True
    return False
e = np.zeros(n, dtype=object); e[0] = 1; add(e); add(R1); add(R2)
i = 0
while i < len(vecs):
    x = vecs[i]; i += 1
    add(gconv(x, R1)); add(gconv(x, R2))
p2 = tuple(int(v[idx[(0,2)]]) for v in vecs); p9 = tuple(int(v[idx[(2,0)]]) for v in vecs)
allok &= ok(len(vecs) == 11 and p2 == p9, 'arc algebra fails to separate: dimT=11, T-profiles of merged pair equal')

# ---- arc-rigidity refutation: the 10-regular thin-merging instance ----
S1 = [(0,1),(0,7),(1,0),(1,2),(2,0),(2,2),(2,4),(2,6),(3,0),(3,6)]
A = cay(S1); g = wl2(A); l = ts_partition(A)
c1 = int(g[0, idx[(2,0)]]); c2 = int(g[0, idx[(2,2)]])
B1 = {els[y] for y in range(n) if g[0, y] == c1}; B2 = {els[y] for y in range(n) if g[0, y] == c2}
thin = B1 <= set(S1) and B2 <= set(S1)
merged = l[0, idx[(2,0)]] == l[0, idx[(2,2)]] and c1 != c2
allok &= ok(g.max()+1 == 13 and l.max()+1 == 12 and thin and merged and fused_cert(A, g, c1, c2),
            'arc rigidity refuted: 10-regular instance merges two THIN classes, certified')

# ---- the instance is schurian: Aut has 13 orbitals = 13 coherent classes ----
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
A0 = cay(S0); G0 = nx.from_numpy_array(A0)
autos = [tuple(m[i] for i in range(n)) for m in GraphMatcher(G0, G0).isomorphisms_iter()]
lab = -np.ones((n, n), dtype=np.int64); c = 0
for x in range(n):
    for y in range(n):
        if lab[x, y] >= 0: continue
        for a in autos: lab[a[x], a[y]] = c
        c += 1
g0 = wl2(A0)
joint = len(set(zip(lab.ravel().tolist(), g0.ravel().tolist())))
allok &= ok(len(autos) == 128 and c == 13 and joint == 13,
            f'instance is schurian: |Aut|={len(autos)}, {c} orbitals = {g0.max()+1} coherent classes (schurian)')

# ---- the nine structure-constant violations, and their location ----
f = g0.copy(); f[f == int(g0[0, idx[(2,0)]])] = int(g0[0, idx[(0,2)]])
_, f = np.unique(f, return_inverse=True); f = f.reshape(n, n)
blocks = [(f == c).astype(np.int64) for c in range(f.max()+1)]
viol = []
for a in range(len(blocks)):
    for b in range(len(blocks)):
        P = blocks[a] @ blocks[b]
        for c2 in range(len(blocks)):
            v = P[blocks[c2] > 0]
            if v.size and (v != v[0]).any(): viol.append((a, b, c2))
merged = {int(f[0, idx[e]]) for e in [(0,2),(0,6),(2,0),(2,4)]}
allok &= ok(len(viol) == 9 and {c2 for _, _, c2 in viol} == merged and
            len({a for a, _, _ in viol}) == 3,
            f'certificate: exactly {len(viol)} violations, all on the merged block, from 3 blocks')

# ---- census spot check via group-level code (validated in schur.validate) ----
coh, ts = closures(dims, S1)
allok &= ok(coh.max()+1 == 13 and ts.max()+1 == 12, 'group-level closures agree on census instance')

# ---- level-<=2 separation across the corpus (sec on the level hierarchy) ----
import networkx as nx
from networkx.generators.atlas import graph_atlas_g
def level2_separates(A):
    g = wl2(A); K = g.max()+1
    N = A.shape[0]
    _, mats, types = classes(A)
    thin = [M for M, t in zip(mats, types) if t == 'thin']
    P = np.eye(N, dtype=object); walks = []
    for _ in range(min(N, 14)):
        walks.append(P.copy()); P = P @ A
    feats = (walks + [T1 @ T2 for T1 in thin for T2 in thin]
             + [M @ A for M in mats] + [A @ M for M in mats])
    reps = {}
    for x in range(N):
        for y in range(N): reps.setdefault(int(g[x, y]), (x, y))
    prof = {c: tuple(int(F[x, y]) for F in feats) for c, (x, y) in reps.items()}
    return len(set(prof.values())) == K
cnt = tested = fails = 0
for G in graph_atlas_g():
    if G.number_of_nodes() < 2 or not nx.is_connected(G): continue
    cnt += 1
    if '--full' not in sys.argv and cnt % 25: continue
    tested += 1
    if not level2_separates(nx.to_numpy_array(G).astype(np.int64)): fails += 1
allok &= ok(fails == 0, f'level-<=2 separation: no failure on {tested} of {cnt} connected atlas graphs'
                        f' ({"full" if "--full" in sys.argv else "sample"})')
for name, A_ in [('K4', cfi(list(nx.complete_graph(4).edges()), 4)[0]),
                 ('K33', cfi(list(nx.complete_bipartite_graph(3, 3).edges()), 6)[0]),
                 ('prism', cfi(list(nx.circular_ladder_graph(3).edges()), 6)[0])]:
    allok &= ok(level2_separates(A_), f'level-<=2 separation on CFI({name}) closure')

print('\nALL CHECKS PASSED' if allok else '\n*** SOME CHECKS FAILED')
sys.exit(0 if allok else 1)
