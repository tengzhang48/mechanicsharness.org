import numpy as np
from mod60_hi import mod60


if __name__ == "__main__":
    Ns = [480, 960, 1440, 1920, 2400]
    vals = []
    for N in Ns:
        v = mod60(N, 'chain')
        vals.append(v)
        print("N=%5d  chain lam~60 = %.7f" % (N, v))
    Ns = np.array(Ns, float)
    vals = np.array(vals)
    A = np.vstack([np.ones(len(Ns)), 1 / Ns, 1 / Ns ** 2]).T
    s, *_ = np.linalg.lstsq(A, vals, rcond=None)
    print("extrapolate L + a/N + b/N^2 : L = %.6f" % s[0])