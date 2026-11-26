"""Compare the closed formula with all 136 published terms of OEIS A400778 (antidiagonals 1..16),
and check the five conjectured row formulas and the conjectured inequality far beyond the published range."""
from chain import formula
from functools import lru_cache
F = lru_cache(None)(formula)
data = [int(x) for x in open("oeis_b400778.txt").read().strip().split(",")]
i = 0; bad = 0
for d in range(1, 17):            # antidiagonal d: (n,k) with n+k = d+1
    for n in range(1, d+1):
        k = d+1-n
        if F(min(n,k), max(n,k)) != data[i]: bad += 1; print("MISMATCH", n, k)
        i += 1
print("OEIS terms compared:", i, "mismatches:", bad)

K = 40
T = lambda n, k: F(min(n,k), max(n,k))
A001353 = [0, 1]
for j in range(2, K+2): A001353.append(4*A001353[-1] - A001353[-2])
c1 = all(T(2,k) == 4**(k-1) for k in range(1, K+1))
c2 = all(T(3,k) == 4**(k-1)*A001353[k] for k in range(1, K+1))
c3 = all(T(4,k) == 64*T(4,k-1) - 512*T(4,k-2) for k in range(3, K+1))
c4 = all(T(5,k) == 256*T(5,k-1) - 11776*T(5,k-2) + 65536*T(5,k-3) - 65536*T(5,k-4) for k in range(5, K+1))
c5 = all(T(6,k) == 1024*T(6,k-1) - 237568*T(6,k-2) + 12582912*T(6,k-3) - 150994944*T(6,k-4) for k in range(5, K+1))
print("row formulas n=2..6 hold up to k =", K, ":", c1, c2, c3, c4, c5)
c6 = True
for n in range(1, 13):
    for k in range(n, 13):
        t = T(n,k); bound = 4**((n-1)*(k-1))
        if t > bound or ((t == bound) != (n <= 2 or k <= 2)): c6 = False
print("inequality with equality iff min(n,k)<=2, for n,k<=12:", c6)
