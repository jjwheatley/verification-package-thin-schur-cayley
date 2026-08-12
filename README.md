# Verification code

Computational companion to *Thin-Schur closures of Cayley graphs: positive
circulant regimes, a minimal abelian counterexample, and its census*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21841858.svg)](https://doi.org/10.5281/zenodo.21841858)

Every computational claim in the paper is reproduced here, in exact
arithmetic: integer matrix products throughout, ranks computed independently
modulo two large primes (99999989 and 99999971) with agreement required, and
cyclotomic arithmetic carried out symbolically. No claim rests on
floating-point comparison.

## Requirements

Python 3.10 or later, numpy, networkx; nothing else.

    pip install -r requirements.txt

## Running it

    python verify.py           # ~4-5 minutes
    python verify.py --full    # adds the complete 995-graph atlas sweep

The script prints a PASS or FAIL line for every claim and exits nonzero if any
fails. A successful run ends with `ALL CHECKS PASSED`.

    python schur.py            # cross-validates the group-level code

## What is checked

- **Walk separation at prime order** (positive theorem): pairwise distinct
  walk profiles across the coherent classes of C13(1,5), C17(1,2,4,8),
  C19(1,7,8).
- **The degenerate spectrum of C25(1,2)**: 13 classes, only 12 distinct
  eigenvalues, walk separation nonetheless.
- **The three wedge witnesses** (n = 24, 40, 48): the identity holds while
  walk data does not separate, with the walk tie, the (A∘A³)A values, and the
  weighted-adjacency weights that recover the arc classes.
- **The Bannai-Muzychuk cardinality gap** on C24(1,5), in exact cyclotomic
  arithmetic: dim T = 7 but |Lambda| = 9, and the dual-side algebra has
  dimension 9, so the candidate fusion is not one.
- **The counterexample** over Z4 x Z8: K = 13, dim TS = 12, and the stall
  certificate verified directly on the fused partition.
- **The nine structure-constant violations**: exactly nine, all on the merged
  block, arising from three blocks.
- **The instance is schurian**: |Aut| = 128 with 13 orbitals, matching the 13
  coherent classes.
- **The arc algebra fails to separate**: dim T = 11, equal profiles on the
  merged pair.
- **Arc rigidity is refuted**: a 10-regular census instance merging two thin
  classes, certified.
- **Level-<=2 separation** across the atlas corpus and the CFI closures over
  K4, K3,3 and the prism.
- **Agreement** between the group-level and matrix-level closures.

## Scan data

The two exhaustive computations are recorded under `outputs/`, so the counts
can be checked without rerunning them. Each JSON records, per group, the
number of connected instances tested and every failure found.

- `outputs/minimality/` — all 44 abelian groups of order below 32, summing to
  259,469 connected instances with no failures.
- `outputs/census/` — Z4 x Z8 alone, summing to 129,600 connected instances
  with 512 counterexamples.

The directories correspond to the paper's two separate claims and each sums to
the figure quoted there.

To reproduce them (about fifteen minutes in total):

    python scan_minimality.py {small|mid|high|cyclic1|cyclic2|cyclic3|g48}

Scan output is written to the current directory rather than to `outputs/`, so
a rerun produces a local file to compare against rather than overwriting the
reference copies. `g48` corresponds to `outputs/census/`, the rest to
`outputs/minimality/`. Cyclic orders 4-16 and the Z4 x Z8 census by mask range
are run through `scan_group` directly, as described in the scan script.

Orbit analysis of the 512 counterexamples under Aut(Z4 x Z8) is a short
post-processing step: enumerate the 128 automorphisms as the maps determined
by generator images, then canonicalise each failing connection set under them.

## Files

- `kwl.py` — coherent closure (2-WL), the neighbour-local refinement, and the
  CFI construction.
- `certificates.py` — class decomposition, stall certificates, merged-pair
  location, and exact two-prime rank computations.
- `schur.py` — group-level computation of both closures for abelian Cayley
  graphs, valid because both partitions are translation-invariant;
  cross-validated against the matrix-level code.
- `scan_minimality.py` — exhaustive scans over all inverse-closed connection
  sets of a given abelian group.
- `verify.py` — the verification script described above.
- `The general counterexamples` — (34, 38, 42 vertices): K = 120, 156, 97 with
  dim TS one less in each case, and the stall certificate verified directly.
## Citing

See `CITATION.cff`, or the "Cite this repository" link on the repository page.

## Licence

MIT.
