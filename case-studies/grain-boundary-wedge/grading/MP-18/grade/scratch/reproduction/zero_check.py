import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast
from fd_nonuniform import build_D1D2


def check_zero_mode(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D3 = build_D3_nonuniform(nodes)
    M = M_phys_fast(x, nodes)
    A = -(M @ D3)
    w, V = np.linalg.eig(A)
    i = np.argmax(w.real)     # the ~0 mode (max real part)
    v = V[:, i].real
    v = v / v[-1]
    # continuum residual of the ODE form: does  M D3 v = 0 ?
    Mapply = M @ (D3 @ v)
    print("N=%d zero-mode rate=%.3e  ||M D3 v||=%.3e  ||v||=%.3f  f(0.5)=%.4f" %
          (N, w[i].real, np.linalg.norm(Mapply), np.linalg.norm(v), np.interp(0.5, x, v)))
    return x, v


if __name__ == "__main__":
    for N in [40, 80, 160, 320]:
        check_zero_mode(N)