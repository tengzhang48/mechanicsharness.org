"""Independent over-determined check: refine the local minima of sigma_min(A-lam M).

A is (N-1) x (N-2), M is (N-1) x (N-2): every interior equation is used.
sigma_min does not reach exactly 0 (the redundant equation is inconsistent at
O(discretisation error)) but its local minima sit at the eigenvalues.
"""
import sys
import numpy as np
from scan2 import pencil, smin, refine

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    x, A, M = pencil(N, ngl=700)
    for a, b in [(9.0, 14.0), (255.0, 275.0)]:
        r, s = refine(A, M, a, b)
        print("N=%d  refined minimum: lam=%.5f  smin=%.3e" % (N, r, s))