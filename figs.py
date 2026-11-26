"""Generate the TikZ pictures of the paper (mirror curves traced from the definition) and check the
spanning-tree form of the main formula:  T(n,k)^2 = 4^g * tau(D),  D = corners joined diagonally, boundary merged."""
import itertools, random
from chain import det_int, formula

OPP = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E'}
def mid(c, d):
    i, j = c
    return {'N': (i-.5, j), 'E': (i, j-.5), 'S': (i-.5, j-1), 'W': (i-1, j-.5)}[d]
def nb(c, d):
    i, j = c
    return {'N': (i, j+1), 'S': (i, j-1), 'E': (i+1, j), 'W': (i-1, j)}[d]

def curves(n, k, mirrors):
    """mirrors: set of frozenset({cell, cell}).  Returns the closed curves as lists of midpoints."""
    cells = {(i, j) for i in range(1, n+1) for j in range(1, k+1)}
    def order(c): return ['N', 'E', 'S', 'W'] if sum(c) % 2 == 0 else ['N', 'W', 'S', 'E']
    def nxt(c, d): o = order(c); return o[(o.index(d)+1) % 4]
    arcs = [(c, d) for c in sorted(cells) for d in 'NESW']        # arc leaving side d of cell c
    seen, out = set(), []
    for a in arcs:
        if a in seen: continue
        path, cur = [], a
        while cur not in seen:
            seen.add(cur); c, d = cur
            path.append(mid(c, d))
            e = nxt(c, d)                                           # arc (c,d) ends at side e of c
            c2 = nb(c, e)
            if c2 in cells and frozenset((c, c2)) not in mirrors: cur = (c2, OPP[e])
            else: cur = (c, e)
        out.append(path)
    return out

def internal(n, k):
    s = []
    for i in range(1, n+1):
        for j in range(1, k+1):
            if i < n: s.append(frozenset(((i, j), (i+1, j))))
            if j < k: s.append(frozenset(((i, j), (i, j+1))))
    return s

def tikz(n, k, mirrors, scale=0.62, colors=('black', 'blue!75!black', 'red!75!black', 'green!50!black')):
    L = [r"\begin{tikzpicture}[scale=%s,line cap=round]" % scale,
         r"\draw[gray!55,thin] (0,0) rectangle (%d,%d);" % (n, k)]
    for s in mirrors:
        (a, b) = sorted(s)
        if a[0] != b[0]:   # vertical side x = a[0], between y=a[1]-1 and a[1]
            L.append(r"\draw[line width=1.6pt,gray!70] (%d,%.2f) -- (%d,%.2f);" % (a[0], a[1]-0.86, a[0], a[1]-0.14))
        else:
            L.append(r"\draw[line width=1.6pt,gray!70] (%.2f,%d) -- (%.2f,%d);" % (a[0]-0.86, a[1], a[0]-0.14, a[1]))
    for t, p in enumerate(curves(n, k, mirrors)):
        pts = " -- ".join("(%.1f,%.1f)" % q for q in p)
        L.append(r"\draw[%s,line width=0.9pt,rounded corners=%.1fpt] %s -- cycle;" % (colors[t % len(colors)], 9*scale, pts))
    for i in range(1, n+1):
        for j in range(1, k+1):
            L.append(r"\fill (%.1f,%.1f) circle (1.3pt);" % (i-.5, j-.5))
    L.append(r"\end{tikzpicture}")
    return "\n".join(L)

def sym_search(n, k, want_mirrors, seed=1):
    """A mirror set invariant under both reflections of the rectangle with a single curve."""
    S = internal(n, k)
    def img(s, f): return frozenset(f(c) for c in s)
    fx = lambda c: (n+1-c[0], c[1]); fy = lambda c: (c[0], k+1-c[1])
    orbits, seen = [], set()
    for s in S:
        if s in seen: continue
        o = {s, img(s, fx), img(s, fy), img(img(s, fx), fy)}
        seen |= o; orbits.append(o)
    rnd = random.Random(seed); best = None
    for _ in range(20000):
        pick = [o for o in orbits if rnd.random() < 0.35]
        M = set().union(*pick) if pick else set()
        if len(M) in want_mirrors and len(curves(n, k, M)) == 1:
            return M
    return best

def check_tree_form(nmax=7):
    for n in range(2, nmax+1):
        for k in range(n, nmax+1):
            pts = [(p, q) for p in range(1, n) for q in range(1, k)]
            idx = {v: t for t, v in enumerate(pts)}
            L = [[0]*len(pts) for _ in pts]
            for (p, q) in pts:
                L[idx[(p, q)]][idx[(p, q)]] = 4
                for dp, dq in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    w = (p+dp, q+dq)
                    if w in idx: L[idx[(p, q)]][idx[w]] -= 1
            assert formula(n, k)**2 == 4**((n-1)*(k-1)) * det_int(L), (n, k)
    return True

def check_parity(max_sides=12, nmax=6, samples=300, seed=0):
    """Remark 8 of the paper: the number of curves of a placement X is congruent to gcd(n,k) + |X| modulo 2.
    Checked for every placement on the grids with at most max_sides internal sides, and for random placements
    on all grids with n, k <= nmax."""
    from math import gcd
    for n in range(1, nmax+1):
        for k in range(n, 2*max_sides):
            S = internal(n, k)
            if len(S) > max_sides: break
            for r in range(len(S)+1):
                for X in itertools.combinations(S, r):
                    if (len(curves(n, k, set(X))) - gcd(n, k) - r) % 2: return False
    rnd = random.Random(seed)
    for n in range(1, nmax+1):
        for k in range(n, nmax+1):
            S = internal(n, k)
            for _ in range(samples):
                X = {x for x in S if rnd.random() < 0.5}
                if (len(curves(n, k, X)) - gcd(n, k) - len(X)) % 2: return False
    return True

if __name__ == "__main__":
    print("spanning-tree form verified:", check_tree_form())
    print("parity remark verified:", check_parity())
    A = set()                                               # 3x3 without mirrors
    S3 = internal(3, 3)
    two = next(set(c) for c in itertools.combinations(S3, 2) if len(curves(3, 3, set(c))) == 1)
    C = sym_search(5, 4, range(8, 13))
    print("3x3 no mirrors:", len(curves(3, 3, A)), "curves; with 2 mirrors:", len(curves(3, 3, two)),
          "; 5x4 symmetric design mirrors:", len(C), "curves:", len(curves(5, 4, C)))
    # count check on 3x3: number of single-curve placements must be 240, of which 48 with two mirrors
    cnt = sum(1 for r in range(len(S3)+1) for c in itertools.combinations(S3, r) if len(curves(3, 3, set(c))) == 1)
    cnt2 = sum(1 for c in itertools.combinations(S3, 2) if len(curves(3, 3, set(c))) == 1)
    print("tracer count T(3,3) =", cnt, "; with exactly two mirrors:", cnt2)
    open("fig_a.tex", "w").write(tikz(3, 3, A))
    open("fig_b.tex", "w").write(tikz(3, 3, two))
    open("fig_c.tex", "w").write(tikz(5, 4, C))
