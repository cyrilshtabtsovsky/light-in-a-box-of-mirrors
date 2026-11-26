"""
Numerical/exact checks of the individual steps of the proof, in the numbering of the paper (mirror_designs.pdf).

 L1  Lemma 9: F' (vertices = all side midpoints, arcs = oriented diamond arcs) is the line digraph of the corner
     digraph Gamma, and every side is exactly one arc of Gamma.
 L2  Lemma 10 and Proposition 11 (exact): for a sample of arcs s0 = (u0 -> v0) of Gamma, t_{s0}(F') = T(n,k) and
     t_{s0}(F') = 2^g * t_{u0}(Gamma) / outdeg(v0).  Also Corollary 4 with the root (1,0).
 L3  A count equivalent to equation (3): a + b = g + n + k - 2, where a = #even corners with 0<q<k and
     b = #odd corners with 0<p<n are the numbers of corners of out-degree 2.
 L4  Lemmas 12 and 13 (floating point): the spectrum of Xi = D_O^{-1} (D_O - C D_E^{-1} B) is
     {1 - cos(l pi/n) cos(m pi/k)}, one per orbit of (l,m)~(n-l,k-m); right null vector 1, left null vector
     dh(p) dv(q); the sum of dh*dv over the odd corners is 2nk.
"""
from fractions import Fraction
from chain import det_int, corner_digraph, formula
import numpy as np, math, sys

def cells_arcs(n, k):
    """Arcs of F': for each cell (i,j), 1<=i<=n, 1<=j<=k, occupying [i-1,i]x[j-1,j], the four diamond arcs,
    oriented clockwise if i+j is even and counter-clockwise if i+j is odd.  A side is named by its two corners."""
    arcs = []
    for i in range(1, n+1):
        for j in range(1, k+1):
            top = ((i-1, j), (i, j)); right = ((i, j-1), (i, j)); bottom = ((i-1, j-1), (i, j-1)); left = ((i-1, j-1), (i-1, j))
            cyc = [top, right, bottom, left] if (i+j) % 2 == 0 else [top, left, bottom, right]
            for t in range(4):
                arcs.append((cyc[t], cyc[(t+1) % 4]))
    return arcs

def check(n, k):
    g = (n-1)*(k-1)
    V, garcs = corner_digraph(n, k)
    # --- L1
    side_of = {}
    for (a, b) in garcs:
        side_of[(a, b)] = tuple(sorted((V[a], V[b])))
    line_arcs = set()
    for (a, b) in garcs:
        for (c, d) in garcs:
            if b == c: line_arcs.add((side_of[(a, b)], side_of[(c, d)]))
    fa = cells_arcs(n, k)
    assert len(fa) == len(set(fa)) == 4*n*k
    assert set(fa) == line_arcs, "F' is not the line digraph of Gamma"
    assert len(set(side_of.values())) == len(garcs) == 2*n*k + n + k      # every side is exactly one arc
    # --- L2 (exact)
    sides = sorted(set(side_of.values())); sidx = {s: i for i, s in enumerate(sides)}
    M = len(sides)
    LF = [[0]*M for _ in range(M)]
    for (s, t) in fa:
        LF[sidx[s]][sidx[s]] += 1; LF[sidx[s]][sidx[t]] -= 1
    N = len(V); out = [0]*N
    LG = [[0]*N for _ in range(N)]
    for (a, b) in garcs:
        out[a] += 1; LG[a][a] += 1; LG[a][b] -= 1
    T = formula(min(n, k), max(n, k))
    for (u0, v0) in garcs[::max(1, len(garcs)//5)]:
        s0 = sidx[side_of[(u0, v0)]]
        tF = det_int([[LF[i][j] for j in range(M) if j != s0] for i in range(M) if i != s0])
        tG = det_int([[LG[i][j] for j in range(N) if j != u0] for i in range(N) if i != u0])
        assert tF == T, (tF, T)
        assert Fraction(2**g * tG, out[v0]) == tF
    # --- Corollary 4: T = 2^g * (number of spanning trees of Gamma oriented towards the corner (1,0))
    u = V.index((1, 0))
    assert 2**g * det_int([[LG[i][j] for j in range(N) if j != u] for i in range(N) if i != u]) == T
    # --- L3
    a = sum(1 for (p, q) in V if (p+q) % 2 == 0 and 0 < q < k)
    b = sum(1 for (p, q) in V if (p+q) % 2 == 1 and 0 < p < n)
    assert a + b == g + n + k - 2
    # --- L4 (floating point)
    E = [v for v in V if sum(v) % 2 == 0]; O = [v for v in V if sum(v) % 2 == 1]
    dh = lambda p: 1 if p in (0, n) else 2
    dv = lambda q: 1 if q in (0, k) else 2
    B = np.array([[1.0 if (e[0] == o[0] and abs(e[1]-o[1]) == 1) else 0.0 for o in O] for e in E])
    C = np.array([[1.0 if (e[1] == o[1] and abs(e[0]-o[0]) == 1) else 0.0 for e in E] for o in O])
    DE = np.diag([float(dv(e[1])) for e in E]); DO = np.diag([float(dh(o[0])) for o in O])
    S = DO - C @ np.linalg.inv(DE) @ B
    Xi = np.linalg.inv(DO) @ S
    ev = np.sort(np.linalg.eigvals(Xi).real)
    pred, seen = [], set()
    for l in range(n+1):
        for m in range(k+1):
            if (l, m) in seen or (l, m) == (n-l, k-m): continue
            seen.add((l, m)); seen.add((n-l, k-m))
            pred.append(1 - math.cos(l*math.pi/n)*math.cos(m*math.pi/k))
    pred = np.sort(np.array(pred))
    assert len(pred) == len(O) and np.allclose(ev, pred, atol=1e-8), "spectrum mismatch"
    ell = np.array([dh(o[0])*dv(o[1]) for o in O], dtype=float)
    assert np.allclose(ell @ Xi, 0, atol=1e-9) and np.allclose(Xi @ np.ones(len(O)), 0, atol=1e-9)
    assert abs(ell.sum() - 2*n*k) < 1e-9
    return True

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in range(1, nmax+1):
        for k in range(1, nmax+1):
            check(n, k)
    print(f"all lemma checks passed for 1 <= n,k <= {nmax}")
