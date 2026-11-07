exec(open('<round-runs>/runs/MP-40/explore38.py').read().split('\nprint("N   n')[0])

import mpmath as mp
mp.mp.dps = 130

print("dps=130")
print("N   n   spectrum(first 6)", flush=True)
for N in [16, 20, 24, 28, 32, 36, 40]:
    basis, nodes, L, ev, cand, d = study2(N)
    print("%2d %3d  %s" % (N, len(basis), "  ".join("%.5f" % c[0] for c in cand[:6])), flush=True)