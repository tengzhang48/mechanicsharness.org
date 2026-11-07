import numpy as np
from conv_phys import leading


def slow_mode(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    lam, ev, M, D1, D2 = leading(nodes)
    A = -(M @ D1 @ D2)
    w, V = np.linalg.eig(A)
    i = np.argmax(w.real)
    v = V[:, i].real
    v = v / v[-1]
    return x, -w[i].real, v


if __name__ == "__main__":
    N = 320
    x, lam, v = slow_mode(N)
    print("lambda1 =", lam)
    # inspect near x=0: does v ~ x ln x  or x  or x^2 ?
    print("x        v            v/x        v/(x ln x)")
    for xv in [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.5, 0.8, 0.95]:
        vv = np.interp(xv, x, v)
        print("%.4f  %+.6f   %+.4f   %+.4f" %
              (xv, vv, vv / xv, vv / (xv * np.log(xv + 1e-300))))