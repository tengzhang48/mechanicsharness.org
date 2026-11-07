import numpy as np
import scipy.linalg as sla
from fd_nonuniform import build_D1D2
from phys_oper import M_phys


def march_decay(nodes, T=20.0, nsteps=4000, seed=None):
    """Implicit march of  u_tau = D2 (1 - M D1 u)  with physical kernel.
    Returns times and ||s||_2, and the fitted late-time rate."""
    xcols = np.asarray(nodes[1:-1], float)
    n = len(xcols)
    D1, D2 = build_D1D2(nodes)
    M = M_phys(xcols, nodes)
    u = np.zeros(n)
    dt = T / nsteps
    A = np.eye(n) - dt * (D2 @ M @ D1)
    rhs_const = -dt * (D2 @ np.ones(n))     # u_tau = D2(1 - M D1 u) = D2 - D2 M D1 u
    lu = sla.lu_factor(A)
    ts = np.zeros(nsteps + 1)
    ns = np.zeros(nsteps + 1)
    def compute_s(u):
        return 1.0 - M @ (D1 @ u)
    ns[0] = np.linalg.norm(compute_s(u))
    for k in range(nsteps):
        u = sla.lu_solve(lu, u + rhs_const)
        ts[k + 1] = (k + 1) * dt
        ns[k + 1] = np.linalg.norm(compute_s(u))
    return ts, ns


if __name__ == "__main__":
    for N in [40, 80, 160]:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        ts, ns = march_decay(nodes, T=20.0, nsteps=4000)
        lo = slice(int(len(ts) * 0.6), len(ts))
        rate = -np.polyfit(ts[lo], np.log(ns[lo]), 1)[0]
        print("N=%3d  late-time decay rate = %.6f  ||s||(0,1,5,20)=%s" %
              (N, rate, np.round([ns[np.argmin(np.abs(ts - t))] for t in [0, 1, 5, 20]], 5)))