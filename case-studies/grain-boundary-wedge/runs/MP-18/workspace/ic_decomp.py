import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast


def ic_decomposition(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    D3 = build_D3_nonuniform(nodes)
    A = -(M @ D3)
    w, V = np.linalg.eig(A)
    # IC: s(x,0)=1 on interior nodes (with s(0)=0 imposed via the operator).
    s0 = np.ones(N)
    coef = np.linalg.solve(V, s0)        # s0 = V @ coef
    r = -w.real
    # amplitude magnitude weighted by decay: at large tau, component ~ coef_k e^{-r_k tau}
    idx = np.argsort(r)                  # slowest first
    print("N=%d  mode decomposition of IC s=1 (slowest 6):" % N)
    for k in idx[:6]:
        print("   rate=%10.4f  |coef|=%.4e  arg coef=%.3f" %
              (r[k], np.abs(coef[k]), np.angle(coef[k])))
    # contribution at tau=0 (should sum to ||s0||)
    print("   sum |coef|^2 = %.4f  ||s0||^2=%d" %
          (np.sum(np.abs(coef) ** 2), N))


if __name__ == "__main__":
    ic_decomposition(120)