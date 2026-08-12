"""Exhaustive scan for TS < Coh over abelian Cayley graphs.

Enumerates every inverse-closed connection set (no size cap) over the given
groups, filters to connected (S generates), and computes both closures at the
group level (schur.closures).  Records every instance with dim TS < K.
"""
import numpy as np, itertools, json, sys, time
from schur import _conv_table, closures

def involutions_pairs(dims, els, neg):
    n = len(els)
    zero = 0
    inv = []; pairs = []
    seen = set()
    for i in range(1, n):
        if i in seen: continue
        j = int(neg[i])
        if j == i: inv.append(i)
        else:
            seen.add(i); seen.add(j); pairs.append((i, j))
    return inv, pairs

def generates(members, add, n):
    seen = {0}
    frontier = list(members)
    for m in frontier:
        if m not in seen: seen.add(m)
    frontier = list(seen)
    while frontier:
        new = []
        for a in frontier:
            for m in members:
                b = int(add[a, m])
                if b not in seen:
                    seen.add(b); new.append(b)
        frontier = new
    return len(seen) == n

def scan_group(dims, out, limit=None, mask_range=None):
    tables = _conv_table(dims)
    els, idx, add, neg = tables
    n = len(els)
    inv, pairs = involutions_pairs(dims, els, neg)
    total = 2 ** (len(inv) + len(pairs))
    t0 = time.time(); tested = 0; fails = []
    lo, hi = (1, total) if mask_range is None else mask_range
    for mask in range(lo, hi):
        members = []
        for b, i in enumerate(inv):
            if mask >> b & 1: members.append(i)
        for b, (i, j) in enumerate(pairs):
            if mask >> (len(inv) + b) & 1: members.extend((i, j))
        if len(members) < 2: continue
        if not generates(members, add, n): continue
        tested += 1
        S = [els[i] for i in members]
        coh, ts = closures(dims, S, add_tables=tables)
        K, T = int(coh.max()) + 1, int(ts.max()) + 1
        if T < K:
            fails.append({'dims': dims, 'S': S, 'K': K, 'dimTS': T})
            print('  FAIL', dims, sorted(S), 'K=%d T=%d' % (K, T), flush=True)
        if limit and tested >= limit: break
    rec = {'dims': dims, 'configs': total - 1, 'connected_tested': tested,
           'failures': fails, 'seconds': round(time.time() - t0, 1)}
    out.append(rec)
    print('group %s: %d connected instances, %d failures (%.0fs)' %
          (dims, tested, len(fails), rec['seconds']), flush=True)

if __name__ == '__main__':
    which = sys.argv[1]
    groups = {
        'small':  [(2,4),(2,2,2),(3,3),(2,6),(2,2,4),(4,4),(2,8),(2,2,2,2)],
        'mid':    [(3,6),(2,10),(2,12),(2,2,6),(5,5)],
        'high':   [(3,9),(3,3,3),(2,14)],
        'cyclic1':[(17,),(18,),(19,),(20,),(21,),(22,),(23,),(24,),(25,)],
        'cyclic2':[(26,),(27,),(28,),(29,),(30,)],
        'cyclic3':[(31,)],
        'g48':    [(4,8)],
        'g216':   [(2,16)],
    }[which]
    out = []
    for dims in groups:
        scan_group(dims, out)
    json.dump(out, open('scan_%s.json' % which, 'w'))
