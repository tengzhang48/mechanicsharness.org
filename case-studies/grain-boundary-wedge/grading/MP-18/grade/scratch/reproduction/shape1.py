import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast


def mode1_shape(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    A = -(M @ build_D3_nonuniform(nodes))
    w, V = np.linalg.eig(A)
    r = -w.real
    r = np.where(r > 1, r, np.inf)
    i = np.argmin(r)
    v = V[:, i].real
    v = v / v[-1]
    return x, v


if __name__ == "__main__":
    print("N       f(0.5)        f(0.25)      f(0.75)     f(0.1)")
    vals = []
    Ns = [120, 240, 480, 960, 1440, 1920]
    for N in Ns:
        x, v = mode1_shape(N)
        f05 = np.interp(0.5, x, v)
        vals.append(f05)
        print("%5d   %.7f   %.7f   %.7f   %.7f" %
              (N, f05, np.interp(0.25, x, v), np.interp(0.75, x, v), np.interp(0.1, x, v)))
    Ns = np.array(Ns, float)
    vals = np.array(vals)
    A = np.vstack([np.ones(3), 1 / Ns[-3:], 1 / Ns[-3:] ** 2]).T
    s, *_ = np.linalg.lstsq(A, vals[-3:], rcond=None)
    print("f(0.5) extrap = %.7f" % s[0])