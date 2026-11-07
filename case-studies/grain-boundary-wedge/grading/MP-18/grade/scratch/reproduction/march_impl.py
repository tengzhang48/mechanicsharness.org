import numpy as np
import scipy.linalg as sla
from fd_nonuniform import build_D1D2
from statement_oper import Smat_statement
import march


def march_u(nodes, Mmat, T=5.0, nsteps=2000, u0=None, D1=None, D2=None):
    """Implicit march of  u_tau = D2 (1 - M D1 u)  on interior nodes.
    Returns times and the residual norm ||s||."""
    xcols = np.asarray(nodes[1:-1], float)
    if D1 is None or D2 is None:
        D1, D2 = build_D1D2(nodes)
    n = len(xcols)
    if u0 is None:
        u0 = np.zeros(n)
    dt = T / nsteps
    A = np.eye(n) - dt * (D2 @ Mmat @ D1)      # for implicit step
    # u_{k+1} = A^{-1}(u_k + dt * D2 @ ones)
    b_const = dt * (D2 @ np.ones(n))
    u = u0.copy()
    ts = [0.0]
    ns = [np.linalg.norm(1.0 - Mmat @ (D1 @ u))]
    for k in range(nsteps):
        u = np.linalg.solve(A, u + b_const)
        ts.append((k + 1) * dt)
        ns.append(np.linalg.norm(1.0 - Mmat @ (D1 @ u)))
    return np.array(ts), np.array(ns)


if __name__ == "__main__":
    N = 60
    xn = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], xn, [1.0]))
    D1, D2 = build_D1D2(nodes)
    # statement kernel
    Ms = Smat_statement(xn, nodes)
    ts, ns = march_u(nodes, Ms, T=5.0, nsteps=2000, D1=D1, D2=D2)
    # fit late-time decay rate
    lo = slice(len(ts) * 3 // 5, len(ts))
    rate = -np.polyfit(ts[lo], np.log(ns[lo]), 1)[0]
    print("STATEMENT kernel: late-time decay rate =", rate)
    print("  ||s|| at tau=0,1,2,5:", [ns[np.argmin(np.abs(ts - t))] for t in [0, 1, 2, 5]])
    # march kernel
    Mm = march.Smat(N)
    ts2, ns2 = march_u(nodes, Mm, T=5.0, nsteps=2000, D1=D1, D2=D2)
    rate2 = -np.polyfit(ts2[lo], np.log(ns2[lo]), 1)[0]
    print("MARCH kernel: late-time decay rate =", rate2)
    print("  ||s|| at tau=0,1,2,5:", [ns2[np.argmin(np.abs(ts2 - t))] for t in [0, 1, 2, 5]])