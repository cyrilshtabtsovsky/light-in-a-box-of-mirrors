"""Run the brute-force simulation (brute.c) on every grid with at most 27 internal sides and compare with the formula."""
import subprocess, sys
from chain import formula
if subprocess.call(["gcc", "-O2", "-o", "brute", "brute.c"]) != 0: sys.exit("could not compile brute.c")
ok = True
for n in range(1, 5):
    for k in range(n, 29):
        M = 2*n*k - n - k
        if M > 27: break
        out = subprocess.check_output(["./brute", str(n), str(k)]).decode().split("=")[1].strip()
        f = formula(n, k)
        good = int(out) == f
        ok = ok and good
        print(f"T({n},{k}): simulation {out}, formula {f}", "OK" if good else "MISMATCH", flush=True)
print("ALL OK" if ok else "FAILURES")
