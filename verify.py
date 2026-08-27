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

# ---- the general counterexamples (34, 38, 42 vertices), certified here ----
from certificates import classes as _classes, cert_check, merged_pair
BASES = {34: ([(0,1),(0,4),(1,2),(1,3),(1,4),(2,3),(3,4)], 5, 120),
         38: ([(0,1),(0,2),(0,5),(1,2),(1,3),(1,4),(3,4),(4,5)], 6, 156),
         42: ([(0,1),(0,2),(0,5),(0,6),(1,4),(2,3),(3,4),(3,6),(4,5)], 7, 97)}
for nexp, (edges, nv, Kexp) in BASES.items():
    Ac, verts = cfi(edges, nv)
    gc = wl2(Ac); lc = ts_partition(Ac); Kc = gc.max()+1
    mp = merged_pair(Ac, verts)
    ids = [i[0] for i in mp[0]]
    cert = cert_check(Ac, [(ids[0], ids[1])])
    cOK = (cert['mult_closed'] and cert['transpose_closed'] and cert['IAJ']
           and cert['coherence_violations'] > 0)
    allok &= ok(bool(Ac.shape[0] == nexp and Kc == Kexp and lc.max()+1 == Kc-1 and cOK),
                f'general counterexample n={nexp}: K={Kc}, dimTS={Kc-1}, stall certified')

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

# ---- appendix claims: blocks self-inverse; stabiliser is Klein four ----
els4 = [(a, b) for a in range(4) for b in range(8)]
APPENDIX_BLOCKS = [
    [(0,0)], [(0,1),(0,3),(0,5),(0,7)], [(0,2),(0,6),(2,0),(2,4)], [(0,4)],
    [(1,0),(1,6),(3,0),(3,2)], [(1,1),(3,7)], [(1,2),(1,4),(3,4),(3,6)],
    [(1,3),(1,7),(3,1),(3,5)], [(1,5),(3,3)], [(2,1),(2,7)], [(2,2),(2,6)],
    [(2,3),(2,5)]]
neg4 = lambda e: ((-e[0]) % 4, (-e[1]) % 8)
allok &= ok(all(sorted(map(neg4, B)) == sorted(B) for B in APPENDIX_BLOCKS),
            'appendix: every block of the fused partition is self-inverse')
Sset0 = set(S0)
stab = []
for a in els4:
    if a[1] % 2: continue
    for b in els4:
        img = {}
        for x, y in els4:
            img[(x, y)] = ((a[0]*x + b[0]*y) % 4, (a[1]*x + b[1]*y) % 8)
        if len(set(img.values())) != n: continue
        if {img[s] for s in S0} == Sset0: stab.append(img)
nontriv = [m for m in stab if any(m[e] != e for e in els4)]
orders = set()
for m in nontriv:
    o, cur = 1, m
    while any(cur[e] != e for e in els4):
        cur = {e: m[cur[e]] for e in els4}; o += 1
    orders.add(o)
allok &= ok(len(stab) == 4 and orders == {2},
            f'rem:schurian: stabiliser has order {len(stab)}, all non-identity elements of order 2 (Klein four)')

# ---- across-layer eigenvalue coincidences (Section 5) ----
# Exhaustive over connected circulants with 4 <= n <= NMAX; the paper's claim
# is for NMAX = 24 (12,043 circulants, 3,946 with a coincidence).  Exact test:
# lambda_a is the coefficient vector of sum x^(as) reduced mod the n-th
# cyclotomic polynomial; equality of lambdas is equality of those vectors.
try:
    import sympy
    from math import gcd as _gcd
    NMAX = 24 if '--full' in sys.argv else 14
    xs = sympy.symbols('x')
    tot = coin = 0
    for m in range(4, NMAX + 1):
        Phi = sympy.Poly(sympy.cyclotomic_poly(m, xs), xs)
        half = list(range(1, m // 2 + 1))
        for r in range(1, len(half) + 1):
            for gens in itertools.combinations(half, r):
                Sc = sorted({g % m for g in gens} | {(-g) % m for g in gens})
                Ac = np.zeros((m, m), dtype=np.int64)
                for i2 in range(m):
                    for sc in Sc: Ac[i2, (i2 + sc) % m] = 1
                if not nx.is_connected(nx.from_numpy_array(Ac)): continue
                tot += 1
                buckets = {}
                for a2 in range(m):
                    v = [0] * m
                    for sc in Sc: v[(a2 * sc) % m] += 1
                    key = tuple(sympy.rem(sympy.Poly(list(reversed(v)), xs), Phi).all_coeffs())
                    buckets.setdefault(key, []).append(a2)
                if any(len({_gcd(a2, m) for a2 in ids}) > 1 for ids in buckets.values()):
                    coin += 1
    if NMAX == 24:
        allok &= ok(tot == 12043 and coin == 3946,
                    f'across-layer coincidences: {coin} of {tot} connected circulants, 4<=n<=24')
    else:
        allok &= ok(tot > 0 and coin > 0,
                    f'across-layer coincidences: {coin} of {tot} connected circulants, 4<=n<={NMAX}'
                    ' (pass --full for the paper\'s 4<=n<=24 figures)')
except ImportError:
    print('SKIP  across-layer coincidence scan (sympy not installed)')

# ---- smallest across-layer instances (Section 5): K_4, then C_6(1,3), C_6(1,2) ----
try:
    import sympy
    from math import gcd as _gcd
    xs2 = sympy.symbols('x')
    def _has_coincidence(m, Sc):
        Phi = sympy.Poly(sympy.cyclotomic_poly(m, xs2), xs2)
        b = {}
        for a2 in range(m):
            v = [0] * m
            for sc in Sc: v[(a2 * sc) % m] += 1
            b.setdefault(tuple(sympy.rem(sympy.Poly(list(reversed(v)), xs2), Phi).all_coeffs()),
                         []).append(a2)
        return any(len({_gcd(a2, m) for a2 in ids}) > 1 for ids in b.values())
    hits = []
    for m in (4, 5, 6):
        half = list(range(1, m // 2 + 1))
        for r in range(1, len(half) + 1):
            for gens in itertools.combinations(half, r):
                Sc = sorted({g % m for g in gens} | {(-g) % m for g in gens})
                Ac = np.zeros((m, m), dtype=np.int64)
                for i2 in range(m):
                    for sc in Sc: Ac[i2, (i2 + sc) % m] = 1
                if not nx.is_connected(nx.from_numpy_array(Ac)): continue
                if _has_coincidence(m, Sc):
                    hits.append((m, tuple(gens), len(Sc) == m - 1))
    smallest = hits[0] if hits else None
    noncomplete = [h for h in hits if not h[2]]
    allok &= ok(smallest == (4, (1, 2), True)
                and sorted(h[1] for h in noncomplete)[:2] == [(1, 2), (1, 3)]
                and all(h[0] == 6 for h in noncomplete),
                'smallest across-layer instance is K_4 = C_4(1,2); smallest '
                'non-complete are C_6(1,3) and C_6(1,2)')
except ImportError:
    print('SKIP  smallest across-layer instances (sympy not installed)')

# ---- cor:lower is vacuous on abelian Cayley graphs (rem:lowerwitness) ----
import random as _rnd
_rnd.seed(0)
bad = 0; seen = 0
for dims in [(4, 8), (16,), (2, 2, 4), (3, 9), (24,)]:
    els2 = list(itertools.product(*[range(d) for d in dims]))
    idx2 = {e: i for i, e in enumerate(els2)}
    nz = [e for e in els2 if any(e)]
    negf = lambda e: tuple((-a) % d for a, d in zip(e, dims))
    for _ in range(8):
        k = _rnd.randint(1, max(1, len(nz) // 2))
        Ss = set()
        for e in _rnd.sample(nz, k): Ss.add(e); Ss.add(negf(e))
        N = len(els2)
        Ac = np.zeros((N, N), dtype=np.int64)
        for i2, xx in enumerate(els2):
            for sc in Ss:
                Ac[i2, idx2[tuple((a + b) % d for a, b, d in zip(xx, sc, dims))]] = 1
        if not nx.is_connected(nx.from_numpy_array(Ac)): continue
        seen += 1
        gg = wl2(Ac)
        if (gg != gg.T).any(): bad += 1
allok &= ok(seen > 0 and bad == 0,
            f'rem:lowerwitness: all coherent classes symmetric on {seen} abelian Cayley instances')

# ---- wedge-degree witnesses: r_3 constant on arcs, r_2 separating ----
for m, gens, r3v, r2v in [(12, [1, 4, 5], 21, {2, 3}), (18, [2, 3, 4, 8], 39, {3, 0})]:
    Ac = np.zeros((m, m), dtype=np.int64)
    for i2 in range(m):
        for gg in gens: Ac[i2, (i2 + gg) % m] = 1; Ac[i2, (i2 - gg) % m] = 1
    gc = wl2(Ac)
    rp = {}
    for y2 in range(m): rp.setdefault(int(gc[0, y2]), y2)
    arcs2 = [c2 for c2, y2 in rp.items() if Ac[0, y2]]
    r3 = Ac * np.linalg.matrix_power(Ac, 3)
    r2 = Ac * np.linalg.matrix_power(Ac, 2)
    allok &= ok(len(arcs2) == 2 and np.array_equal(r3, r3v * Ac)
                and {int(r2[0, rp[c2]]) for c2 in arcs2} == r2v,
                f'rem:wedgedegree C{m}{tuple(gens)}: r_3 = {r3v}A on arcs, r_2 separates')

# ---- prop:strict minimal witness C_20(2,4,5): 3 arc classes, r_2 separates, r_3 does not ----
m = 20; Sw = [2, 4, 5, 15, 16, 18]
Aw = np.zeros((m, m), dtype=np.int64)
for i2 in range(m):
    for sc in Sw: Aw[i2, (i2 + sc) % m] = 1
gw = wl2(Aw); lw = ts_partition(Aw)
rw = {}
for y2 in range(m): rw.setdefault(int(gw[0, y2]), y2)
arcw = [c2 for c2, y2 in rw.items() if Aw[0, y2]]
r2w = [int((Aw * np.linalg.matrix_power(Aw, 2))[0, rw[c2]]) for c2 in arcw]
r3w = [int((Aw * np.linalg.matrix_power(Aw, 3))[0, rw[c2]]) for c2 in arcw]
allok &= ok(len(arcw) == 3 and sorted(r2w) == [0, 1, 2] and sorted(r3w) == [15, 15, 16]
            and gw.max() == lw.max(),
            f'prop:strict minimal witness C_20(2,4,5): 3 arc classes, r_2={sorted(r2w)} separates, '
            f'r_3={sorted(r3w)} does not, AS=Coh')

# ---- prop:wa generalised: Vandermonde recovery with t classes ----
Mw = Aw * np.linalg.matrix_power(Aw, 2)
wts = [int(Mw[0, rw[c2]]) for c2 in arcw]
V = np.array([[1] * len(arcw)] + [[w ** j for w in wts] for j in range(1, len(arcw))], dtype=float)
allok &= ok(len(set(wts)) == len(wts) and abs(np.linalg.det(V)) > 1e-9,
            f'prop:wa: weights {sorted(wts)} pairwise distinct, Vandermonde system inverts (t={len(arcw)})')

# ---- rem:degencyc minimal degenerate cyclotomic witness C_10(1,2) ----
m = 10; Sd = sorted({g % m for g in (1, 2)} | {(-g) % m for g in (1, 2)})
Ad = np.zeros((m, m), dtype=np.int64)
for i2 in range(m):
    for sc in Sd: Ad[i2, (i2 + sc) % m] = 1
gd = wl2(Ad); Kd = gd.max() + 1
zz = np.exp(2j * np.pi / m)
evd = [complex(sum(zz ** ((a2 * sc) % m) for sc in Sd)) for a2 in range(m)]
rounded = [complex(round(e.real, 7), round(e.imag, 7)) for e in evd]
allok &= ok(Kd == 6 and len(set(rounded)) < m,
            f'rem:degencyc: C_10(1,2) has K={Kd} with a repeated orbit eigenvalue (degenerate)')

# ---- rem:counterscan: Z_4 x Z_4 has no counterexample ----
try:
    from schur import closures as _cl, _conv_table as _ct
    tabs44 = _ct((4, 4))
    els44 = tabs44[0]
    neg44 = tabs44[3]
    inv44 = [i2 for i2 in range(1, 16) if int(neg44[i2]) == i2]
    prs44 = []
    seen44 = set()
    for i2 in range(1, 16):
        if i2 in seen44: continue
        j2 = int(neg44[i2])
        if j2 != i2: seen44.update({i2, j2}); prs44.append((i2, j2))
    n44 = tested44 = fails44 = 0
    for mask in range(1, 2 ** (len(inv44) + len(prs44))):
        mem = [inv44[b] for b in range(len(inv44)) if mask >> b & 1]
        for b, (i2, j2) in enumerate(prs44):
            if mask >> (len(inv44) + b) & 1: mem += [i2, j2]
        if len(mem) < 2: continue
        seenset = {0}; frontier = list(mem)
        while frontier:
            nxt = []
            for a2 in frontier:
                for mm in mem:
                    bb = int(tabs44[2][a2, mm])
                    if bb not in seenset: seenset.add(bb); nxt.append(bb)
            frontier = nxt
        if len(seenset) != 16: continue
        tested44 += 1
        co, ts44 = _cl((4, 4), [els44[i2] for i2 in mem], add_tables=tabs44)
        if ts44.max() < co.max(): fails44 += 1
    allok &= ok(tested44 == 432 and fails44 == 0,
                f'rem:counterscan: Z_4 x Z_4 has {tested44} connected instances, {fails44} counterexamples')
except Exception as e:
    print(f'*** FAIL  Z_4 x Z_4 scan raised {e}'); allok = False

# ---- lem:orient parenthetical: AJ at level 1 on irregular closures ----
def _ajsep(Am):
    nn=Am.shape[0]; c=wl2(Am)
    if (c==c.T).all(): return None
    Jm=np.ones((nn,nn),dtype=np.int64); AJ=Am@Jm; JA=Jm@Am
    prs=[(x,y) for x in range(nn) for y in range(nn) if c[x,y]!=c[y,x]]
    good=sum(1 for x,y in prs if AJ[x,y]!=AJ[y,x] or JA[x,y]!=JA[y,x])
    return good,len(prs)
P3=np.array([[0,1,1],[1,0,0],[1,0,0]],dtype=np.int64)
r_p3=_ajsep(P3)
r_k4=_ajsep(cfi(list(nx.complete_graph(4).edges()),4)[0])
r_k33=_ajsep(cfi(list(nx.complete_bipartite_graph(3,3).edges()),6)[0])
r_pr=_ajsep(cfi(list(nx.circular_ladder_graph(3).edges()),6)[0])
allok &= ok(r_p3[0]==r_p3[1] and r_k4[0]==r_k4[1] and r_k33[0]==r_k33[1]
            and r_pr==(864,1008),
            f'lem:orient: AJ separates all transpose pairs on P3, CFI(K4), CFI(K33); '
            f'{r_pr[0]} of {r_pr[1]} on CFI(prism)')

# ---- alternation depth of prime cycles (walk separation is not level 1) ----
def _delta_cycle(m):
    from schur import _conv_table as _ct
    tb=_ct((m,)); el=tb[0]; ix={e:i for i,e in enumerate(el)}; ad=tb[2]
    Sc={(1,), ((-1) % m,)}
    lb=np.full(m,2,dtype=np.int64)
    for e in Sc: lb[ix[e]]=1
    lb[0]=0; dd=0
    while True:
        cc={}
        for i2 in range(m): cc.setdefault(int(lb[i2]),[]).append(i2)
        ft=[lb.copy()]
        for kk,mm in cc.items():
            if 0 not in mm and all(el[i2] in Sc for i2 in mm):
                for kk2,mm2 in cc.items():
                    vv=np.zeros(m,dtype=np.int64)
                    for a2 in mm:
                        for b2 in mm2: vv[int(ad[a2,b2])]+=1
                    ft.append(vv)
        _,nw=np.unique(np.stack(ft,1),axis=0,return_inverse=True)
        if len(set(nw.tolist()))==len(cc): return dd
        lb=nw; dd+=1
_ds={p:_delta_cycle(p) for p in (11,13,17,19,23)}
allok &= ok(_ds=={11:3,13:4,17:6,19:7,23:9},
            f'alternation depth of prime cycles {_ds}: walk separation is not a level-1 phenomenon')

# ---- Bannai-Muzychuk: Sch(V) matrix-closed while |Lambda| != |Delta| ----
_A24=np.zeros((24,24),dtype=np.int64)
for i2 in range(24):
    for sc in (1,5,19,23): _A24[i2,(i2+sc)%24]=1
_g24=wl2(_A24)
from certificates import _rref_add as _rr, P1 as _PP
def _grpalg(gen,mm):
    bas=[]; vs=[]
    def _ad(v):
        vv=np.array([int(x)%_PP for x in v],dtype=np.int64)
        if _rr(bas,vv,_PP): vs.append(v); return True
        return False
    e0=np.zeros(mm,dtype=object); e0[0]=1; _ad(e0)
    for G in gen: _ad(G)
    i3=0
    while i3<len(vs):
        x=vs[i3]; i3+=1
        for G in gen:
            w=np.zeros(mm,dtype=object)
            for a3 in range(mm):
                if x[a3]:
                    for b3 in range(mm):
                        if G[b3]: w[(a3+b3)%mm]+=x[a3]*G[b3]
            _ad(w)
    return vs
_Av=np.array([1 if y2 in {1,5,19,23} else 0 for y2 in range(24)],dtype=object)
_vs=_grpalg([_Av],24)
_pf={}
for y2 in range(24): _pf.setdefault(tuple(int(v[y2]) for v in _vs),[]).append(y2)
allok &= ok(len(_vs)==7 and len(_pf)==int(_g24.max())+1==9,
            f'C24(1,5): dim V={len(_vs)}, dim Sch(V)={len(_pf)}=K, so Sch(V) is matrix-closed '
            'while the BM cardinality condition fails')

# ---- prop:strictinc: Sch(L) is strictly smaller than AS on the CFI witnesses ----
def _dims(Am):
    nn=Am.shape[0]; gg=wl2(Am); ll=ts_partition(Am)
    KK=int(gg.max())+1; dA=int(ll.max())+1
    AI=((Am+np.eye(nn,dtype=np.int64))>0)
    loc=[]
    for c2 in range(dA):
        M=(ll==c2).astype(np.int64)
        if ((M!=0)<=AI).all() and M.any(): loc.append(M)
    def _alg(p):
        bas=[]; vs=[]
        def _ad(M):
            v=(M.reshape(-1)%p).astype(np.int64)
            if _rref_add(bas,v,p): vs.append(M); return True
            return False
        _ad(np.eye(nn,dtype=np.int64))
        for G in loc: _ad(G)
        i3=0
        while i3<len(vs):
            X=vs[i3]; i3+=1
            for G in loc: _ad((X@G)%p); _ad((G@X)%p)
        return vs
    V1=_alg(P1); V2=_alg(P2)
    assert len(V1)==len(V2)
    pr={}
    for x2 in range(nn):
        for y2 in range(nn):
            pr.setdefault(tuple(int(M[x2,y2])%P1 for M in V1),[]).append((x2,y2))
    return KK,dA,len(pr),len(V1)
_exp={34:(120,119,117,114),38:(156,155,151,148),42:(97,96,94,92)}
_got={}
for _n,_bt in BASES.items():
    E2,nv2 = _bt[0],_bt[1]
    Ax,_=cfi(E2,nv2); _got[_n]=_dims(Ax)
Qx,_=cfi(list(nx.convert_node_labels_to_integers(nx.hypercube_graph(3)).edges()),8)
_q=_dims(Qx)
allok &= ok(_got==_exp and _q==(20,20,20,14),
            f'prop:strictinc: (K,dimAS,dimSch(L),dimL) = {_got}, Q_3 {_q}; inclusion strict on all three')

# ---- deposited prop:general certificates are self-consistent ----
import json as _json, os as _os
_cd=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),'outputs','certificates')
_certok=True; _summ=[]
for _nn in (34,38,42):
    _p=_os.path.join(_cd,'cfi_%d.json'%_nn)
    if not _os.path.exists(_p): _certok=False; _summ.append('%d:missing'%_nn); continue
    _r=_json.load(open(_p))
    Ax,_=cfi(_r['base_edges'],_r['base_vertices'])
    gx=wl2(Ax); lx=ts_partition(Ax)
    good=(Ax.shape[0]==_r['n'] and int(gx.max())+1==_r['K']
          and int(lx.max())+1==_r['dim_AS']
          and _r['licensed_products_block_constant'] and _r['transpose_closed']
          and _r['contains_I_A_J'] and _r['structure_constant_violations']>0
          and _r['fused_blocks']==_r['K']-1)
    _certok &= good
    _summ.append('%d:%s'%(_nn,'ok' if good else 'BAD'))
allok &= ok(_certok, 'prop:general certificates deposited and consistent (%s)'%', '.join(_summ))

# ---- rem:cycscope: how much of the prime-power range thm:cyc covers ----
from math import gcd as _gcd2
def _scope(mm):
    zz=np.exp(2j*np.pi/mm); tot=cy=cn=0; degen=[]
    half=list(range(1,mm//2+1))
    for r2 in range(1,len(half)+1):
        for gg in itertools.combinations(half,r2):
            Sx=sorted({g2%mm for g2 in gg}|{(-g2)%mm for g2 in gg})
            Ax=np.zeros((mm,mm),dtype=np.int64)
            for i2 in range(mm):
                for sc in Sx: Ax[i2,(i2+sc)%mm]=1
            if not nx.is_connected(nx.from_numpy_array(Ax)): continue
            Hx={t for t in range(1,mm) if _gcd2(t,mm)==1 and {(t*sc)%mm for sc in Sx}==set(Sx)}
            sn=set(); ob=[]
            for x2 in range(1,mm):
                if x2 in sn: continue
                o2=sorted({(x2*h2)%mm for h2 in Hx}); ob.append(o2); sn|=set(o2)
            Kx=int(wl2(Ax).max())+1
            tot+=1
            if Kx==len(ob)+1:
                cy+=1
                lm=[complex(round(sum(zz**((o2[0]*sc)%mm) for sc in Sx).real,7),
                            round(sum(zz**((o2[0]*sc)%mm) for sc in Sx).imag,7)) for o2 in ob]
                if len(set(lm))==len(lm): cn+=1
                else: degen.append(tuple(Sx))
    return tot,cy,cn,degen
_T=_C=_N=0; _dg=[]
for _m in (4,8,9,16,25,27):
    a,b,c,d=_scope(_m); _T+=a; _C+=b; _N+=c
    if _m==16: _dg=d
allok &= ok(_T==12536 and _C==12454 and _N==10570 and (1,3,13,15) in _dg,
            f'rem:cycscope: {_T} connected prime-power circulants, {_C} cyclotomic, '
            f'{_N} also non-degenerate; C_16(1,3,13,15) is cyclotomic but degenerate')

print('\nALL CHECKS PASSED' if allok else '\n*** SOME CHECKS FAILED')
sys.exit(0 if allok else 1)
