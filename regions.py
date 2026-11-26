"""
Single-line mirror designs on an arbitrary finite set of cells (dots).

Proposition 7 of the paper (mirror placements with a single curve = Eulerian circuits of the oriented arc digraph,
hence = a cofactor of its Laplacian by the BEST theorem) does not use the rectangular shape.
Here it is evaluated for the diamond-shaped dot arrays 1-3-...-(2r-1)-...-3-1 of the kolam literature and
compared with the exhaustive counts 240 and 11,661,312 of Yanagisawa and Nagata (Forma 22 (2007), 31-46).
"""
from chain import det_int
import sys

def euler_count(cells):
    cells = set(cells)
    sides = {}                                   # internal side -> index
    def side(c, d):
        i, j = c
        nb = {'N': (i, j+1), 'S': (i, j-1), 'E': (i+1, j), 'W': (i-1, j)}[d]
        if nb not in cells: return None
        return sides.setdefault(frozenset((c, nb)), len(sides))
    arcs = []
    for c in sorted(cells):
        order = ['N', 'E', 'S', 'W'] if sum(c) % 2 == 0 else ['N', 'W', 'S', 'E']
        cyc = [s for s in (side(c, d) for d in order) if s is not None]
        for t in range(len(cyc)):
            arcs.append((cyc[t], cyc[(t+1) % len(cyc)]))
    m = len(sides)
    if m <= 1: return 1, m
    L = [[0]*m for _ in range(m)]
    for a, b in arcs:
        if a != b: L[a][a] += 1; L[a][b] -= 1
    return det_int([r[1:] for r in L[1:]]), m

def brute(cells):
    """Direct simulation of the definition on an arbitrary cell set (union-find on arcs)."""
    cells = sorted(set(cells)); cs = set(cells); idx = {c: t for t, c in enumerate(cells)}
    internal = []
    for (i, j) in cells:
        if (i+1, j) in cs: internal.append(('v', (i, j), (i+1, j)))
        if (i, j+1) in cs: internal.append(('h', (i, j), (i, j+1)))
    m = len(internal); na = 4*len(cells)       # arcs of a cell: 0=NE 1=ES 2=SW 3=WN
    count = 0
    for x in range(1 << m):
        par = list(range(na))
        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        def uni(a, b):
            a, b = find(a), find(b)
            if a != b: par[a] = b
        for (i, j) in cells:
            A = 4*idx[(i, j)]
            if (i, j+1) not in cs: uni(A+3, A+0)
            if (i, j-1) not in cs: uni(A+1, A+2)
            if (i-1, j) not in cs: uni(A+2, A+3)
            if (i+1, j) not in cs: uni(A+0, A+1)
        for e, (kind, c1, c2) in enumerate(internal):
            A, B = 4*idx[c1], 4*idx[c2]
            mirror = (x >> e) & 1
            if kind == 'v':                      # c1 left, c2 right
                if mirror: uni(A+0, A+1); uni(B+2, B+3)
                else:      uni(A+0, B+2); uni(A+1, B+3)
            else:                                # c1 below, c2 above: c1's N side
                if mirror: uni(A+3, A+0); uni(B+1, B+2)
                else:      uni(A+3, B+1); uni(A+0, B+2)
        r = find(0)
        if all(find(a) == r for a in range(na)): count += 1
    return count, m

def diamond(r):
    return [(i, j) for i in range(-r, r+1) for j in range(-r, r+1) if abs(i) + abs(j) <= r-1]

if __name__ == "__main__":
    for r in range(1, 8):
        t, m = euler_count(diamond(r))
        line = f"diamond 1-...-{2*r-1}-...-1: {len(diamond(r))} dots, {m} internal sides, single-line designs = {t}"
        if m <= 16:
            b, _ = brute(diamond(r)); line += f"   [simulation: {b}]"; assert b == t
        print(line, flush=True)
    # a few non-rectangular shapes, simulation versus determinant
    L_shape = [(0,0),(1,0),(2,0),(0,1),(0,2),(1,1)]
    T_shape = [(0,0),(1,0),(2,0),(1,1),(1,2),(0,2),(2,2)]
    hole = [(i,j) for i in range(3) for j in range(3) if (i,j) != (1,1)]
    for name, cells in [("L", L_shape), ("T/I", T_shape), ("ring", hole)]:
        assert euler_count(cells)[0] == brute(cells)[0], name
        print(name, "shape: determinant = simulation =", euler_count(cells)[0])
