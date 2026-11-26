"""
Verify every step of the proof chain for T(n,k) (OEIS A400778) with exact integer arithmetic.

  (A) T = number of Eulerian circuits of the oriented mirror-curve graph F   [= any cofactor of 2I - A_F, BEST theorem]
  (B) T = 2^g * kappa(Gamma, u0) / outdeg(v0)   for any arc u0->v0 of the corner digraph Gamma
  (C) T^2 = prod_{l=1..n-1} prod_{m=1..k-1} (16 - 16 cos(l pi/n) cos(m pi/k)) = det(16 I - 4 A(P_{n-1}) (x) A(P_{k-1}))
"""
from fractions import Fraction
from math import isqrt
import sys

def det_int(M):
    """Exact determinant of an integer matrix (fraction-free Bareiss)."""
    n = len(M)
    if n == 0: return 1
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for i in range(n-1):
        if M[i][i] == 0:
            for j in range(i+1, n):
                if M[j][i] != 0:
                    M[i], M[j] = M[j], M[i]; sign = -sign; break
            else:
                return 0
        for j in range(i+1, n):
            for c in range(i+1, n):
                M[j][c] = (M[j][c]*M[i][i] - M[j][i]*M[i][c]) // prev
        prev = M[i][i]
    return sign * M[n-1][n-1]

def path_adj(s):
    return [[1 if abs(i-j) == 1 else 0 for j in range(s)] for i in range(s)]

def formula(n, k):
    """T(n,k) from the closed formula, computed exactly as sqrt(det(16 I - 4 A(P_{n-1}) (x) A(P_{k-1})))."""
    a, b = n-1, k-1
    if a == 0 or b == 0: return 1
    A, B = path_adj(a), path_adj(b)
    N = a*b
    M = [[0]*N for _ in range(N)]
    for i in range(a):
        for j in range(b):
            for i2 in range(a):
                for j2 in range(b):
                    M[i*b+j][i2*b+j2] = (16 if (i == i2 and j == j2) else 0) - 4*A[i][i2]*B[j][j2]
    d = det_int(M)
    r = isqrt(d)
    assert r*r == d, "determinant is not a perfect square"
    return r

# ---------- step (A): the oriented mirror-curve graph F ----------
def mirror_digraph(n, k):
    """Vertices: internal sides. Dots (r,c), rows r downward. Each dot's diamond N->E->S->W is a directed
    cycle, clockwise if r+c even and counter-clockwise if r+c odd; border sides are skipped."""
    idx = {}
    def side(r, c, d):   # canonical name of side d in {N,E,S,W} of dot (r,c); None if on the border
        if d == 'N': return None if r == 0 else ('h', r-1, c)
        if d == 'S': return None if r == n-1 else ('h', r, c)
        if d == 'W': return None if c == 0 else ('v', r, c-1)
        if d == 'E': return None if c == k-1 else ('v', r, c)
    sides = [('v', r, c) for r in range(n) for c in range(k-1)] + [('h', r, c) for r in range(n-1) for c in range(k)]
    for i, s in enumerate(sides): idx[s] = i
    arcs = []
    for r in range(n):
        for c in range(k):
            order = ['N', 'E', 'S', 'W'] if (r+c) % 2 == 0 else ['N', 'W', 'S', 'E']
            cyc = [side(r, c, d) for d in order]
            cyc = [s for s in cyc if s is not None]
            if len(cyc) >= 2:
                for i in range(len(cyc)):
                    arcs.append((idx[cyc[i]], idx[cyc[(i+1) % len(cyc)]]))
            # a dot with a single internal side (n=1 or k=1 ends) contributes a loop, which does not affect cofactors
    return len(sides), arcs

def euler_circuits(n, k):
    m, arcs = mirror_digraph(n, k)
    if m <= 1: return 1
    out = [0]*m; ind = [0]*m
    L = [[0]*m for _ in range(m)]
    for a, b in arcs:
        if a == b: continue
        out[a] += 1; ind[b] += 1
        L[a][a] += 1; L[a][b] -= 1
    assert out == ind, "not Eulerian"
    assert all(o <= 2 for o in out)
    return det_int([row[1:] for row in L[1:]])   # BEST: ec = #arborescences * prod (outdeg-1)! = cofactor

# ---------- step (B): the corner digraph Gamma ----------
def corner_digraph(n, k):
    """Corners (p,q), 0<=p<=n (horizontal), 0<=q<=k (vertical). Even corners (p+q even) send an arc along every
    vertical side at them, odd corners along every horizontal side."""
    V = [(p, q) for p in range(n+1) for q in range(k+1)]
    idx = {v: i for i, v in enumerate(V)}
    arcs = []
    for (p, q) in V:
        if (p+q) % 2 == 0: nb = [(p, q-1), (p, q+1)]
        else:              nb = [(p-1, q), (p+1, q)]
        for w in nb:
            if w in idx: arcs.append((idx[(p, q)], idx[w]))
    return V, arcs

def kappa_gamma(n, k):
    V, arcs = corner_digraph(n, k)
    N = len(V)
    out = [0]*N
    L = [[0]*N for _ in range(N)]
    for a, b in arcs:
        out[a] += 1; L[a][a] += 1; L[a][b] -= 1
    vals = set()
    g = (n-1)*(k-1)
    res = []
    for (u, v) in arcs[:6] + arcs[-6:]:
        minor = [[L[i][j] for j in range(N) if j != u] for i in range(N) if i != u]
        kap = det_int(minor)            # arborescences of Gamma directed towards root u
        val = Fraction(2**g * kap, out[v])
        vals.add(val)
    assert len(vals) == 1, vals
    val = vals.pop()
    assert val.denominator == 1
    return int(val)

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    brute = {(1,5):1,(2,2):4,(2,3):16,(2,6):1024,(3,3):240,(3,4):3584,(3,5):53504,(4,4):196608,(2,9):65536,(3,6):798720}
    ok = True
    for n in range(1, nmax+1):
        for k in range(n, nmax+1):
            f = formula(n, k)
            a = euler_circuits(n, k)
            b = kappa_gamma(n, k)
            flag = (f == a == b)
            for key in [(n, k), (k, n)]:
                if key in brute: flag = flag and brute[key] == f
            ok = ok and flag
            print(n, k, f, "OK" if flag else f"MISMATCH euler={a} gamma={b}")
    print("ALL OK" if ok else "FAILURES")
