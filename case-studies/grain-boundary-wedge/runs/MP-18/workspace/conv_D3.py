import numpy as np
from d3_uniform import build_D3_uniform
from phys_oper import M_phys
import march


def leading_D3(N, Mfun=None):
    h = 1.0 / (N + 1)
    x = np.arange(1, N + 1) * h
    nodes = np.concatenate(([0.0], x, [1.0]))
    D3 = build_D3_uniform(N, h)
    if Mfun is None:
        M = M_phys(x, nodes)
    else:
        M = Mfun(x, nodes)
    A = -(M @ D3)
    ev = np.linalg.eigvals(A)
    return x, ev, A, M, D3


if __name__ == "__main__":
    print("=== operator -M D3, physical kernel, uniform grid ===")
    for N in [10, 20, 40, 80, 160]:
        x, ev, A, M, D3 = leading_D3(N)
        print("N=%3d  slowest decay=%.8f  max re=%+.4e  top5 decay=%s" %
              (N, -np.max(ev.real), np.max(ev.real),
               np.round(np.sort(-ev.real)[::-1][:5], 5)))