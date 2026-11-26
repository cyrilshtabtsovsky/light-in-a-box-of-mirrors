import numpy as np, math
G = 0.915965594177219015
lim = math.log(2) + 2*G/math.pi
def logT(n, k):
    l = np.arange(1, n)[:, None]; m = np.arange(1, k)[None, :]
    return (n-1)*(k-1)*math.log(4) + 0.5*np.sum(np.log1p(-np.cos(l*np.pi/n)*np.cos(m*np.pi/k)))
for n in (10, 100, 1000, 4000):
    print(n, logT(n, n)/n**2, "limit", lim, "exp(limit) =", math.exp(lim))
