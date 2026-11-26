# Light in a box of mirrors: counting the placements that give a single closed curve

Proof of the conjectures in OEIS A400778 (single-line mirror designs on an n x k grid).

Cyril Shtabtsovsky (c@astrolabe.vc), 26 November 2026. Worked out with Claude Fable 5.2 (Anthropic).

DOI: [10.5281/zenodo.24546586](https://doi.org/10.5281/zenodo.24546586) (paper and programs, archived on Zenodo)

To cite: C. Shtabtsovsky, *Light in a box of mirrors: counting the placements that give a single closed curve*,
preprint, 2026. https://doi.org/10.5281/zenodo.24546586

## Result

T(n,k) = number of ways to place two-sided mirrors between neighbouring dots of an n x k grid of dots
so that the mirror curve is one closed line. Proved in `mirror_designs.pdf`:

    T(n,k) = 4^((n-1)(k-1)) * prod_{1 <= l < n/2} prod_{1 <= m < k/2} (1 - cos^2(pi*l/n) * cos^2(pi*m/k))

    T(n,k)^2 = prod_{l=1..n-1} prod_{m=1..k-1} (16 - 16*cos(pi*l/n)*cos(pi*m/k))
             = det(16*I - 4*A(P_{n-1}) (x) A(P_{k-1}))        (A(P_r) = adjacency matrix of the path on r vertices)

This implies all six conjectures in the OEIS entry (as it stood on 26 Nov 2026):

1. T(2,k) = 4^(k-1)
2. T(3,k) = 4^(k-1) * A001353(k)
3. T(4,k) = 64 T(4,k-1) - 512 T(4,k-2)                                  (for all k >= 3)
4. T(5,k) = 256 T(5,k-1) - 11776 T(5,k-2) + 65536 T(5,k-3) - 65536 T(5,k-4)      (for all k >= 5)
5. T(6,k) = 1024 T(6,k-1) - 237568 T(6,k-2) + 12582912 T(6,k-3) - 150994944 T(6,k-4)   (for all k >= 5)
6. T(n,k) <= 4^((n-1)(k-1)), with equality exactly when n <= 2 or k <= 2

Also: 2^((n-1)(k-1)) divides T(n,k), and log T(n,k) / (nk) -> log 2 + 2G/pi = 1.27627 (G = Catalan's constant).

## Idea of the proof

1. Orient the diamonds around the dots alternately clockwise / counter-clockwise. A light ray then always
   follows the arrows, with or without a mirror. Single-line designs = Eulerian circuits of this directed graph
   (a piece of the "Manhattan lattice").
2. BEST theorem: number of Eulerian circuits = a determinant (number of oriented spanning trees).
3. The directed graph is the line digraph of a simpler digraph on the grid corners; a theorem of Knuth (1967)
   reduces the determinant to a tree count there.
4. That count is evaluated in closed form with cosine eigenvectors.

## Files

| File | What it is |
|---|---|
| `mirror_designs.pdf` / `.tex` | The paper with the full proof (11 pages) |
| `brute.c` | Direct simulation of the definition, no theory. `gcc -O2 -o brute brute.c && ./brute 3 4` |
| `brute_all.py` | Runs the simulation on every grid with <= 27 internal sides and compares with the formula (about 3 minutes) |
| `chain.py` | Exact check that Euler-circuit count = corner-digraph tree count = formula. `python3 chain.py 8` |
| `lemmas.py` | Checks the steps of the proof separately (needs numpy). `python3 lemmas.py 6` |
| `check_oeis.py` | Compares the formula with all 136 published OEIS terms; tests the six conjectures |
| `asym.py` | Numerical check of the growth constant (needs numpy) |
| `regions.py` | The same count for non-rectangular shapes; reproduces the 2007 kolam counts 240 and 11,661,312 |
| `figs.py` | Traces mirror curves from the definition (used for the mirror-curve pictures); checks the parity remark and the spanning-tree form |
| `oeis_b400778.txt` | The 136 terms of the OEIS b-file, downloaded 25 Nov 2026 |
| `LICENSE`, `CITATION.cff` | MIT for the code, CC BY 4.0 for the paper |

## What was checked by computer (25-26 Nov 2026)

- Simulation = formula for all 41 grids with at most 27 internal sides (largest: 2x9, 3x6, 4x4).
- Formula = all 136 OEIS terms (n + k <= 17).
- Euler-circuit determinant = corner-digraph tree count = formula, exactly, for all n <= k <= 8.
- The steps of the proof (Lemmas 9, 10, 12, 13, Proposition 11), for all n, k <= 6.
- Conjectures 1-5 for k <= 40; conjecture 6 for n, k <= 12.
- The determinant method reproduces the published exhaustive counts for diamond-shaped kolam arrays
  (Yanagisawa-Nagata 2007): 240 and 11,661,312.

## Status

- The proof has not been refereed and has not been formalised in a proof assistant. It is supported by the
  computations above.
- Use of AI: the problem was chosen by the author; the formula, the proof, the programs and the text were
  produced with Claude (Anthropic), model Fable 5.2. See the statement at the end of the paper.
- Prior work (searched 25-26 Nov 2026): earlier work counted single-loop designs only by exhaustive search
  (Yanagisawa-Nagata 2007, diamond arrays); no formula for rectangles was found in the OEIS, arXiv or the web.
  Kasteleyn (1963) and Duplantier-David (1988) count Hamiltonian walks on the Manhattan lattice, a related
  problem; their full texts were not available to us, so the paper says "to our knowledge".
- The 18 references of the paper were checked against publisher, arXiv or index pages on 25-26 Nov 2026.
