import numpy as np
from mfast import M_phys_fast
import scipy.linalg as sla


def transform_op(N, p):
    """Discretise s_tau = -M D3 s  using x = eta^p on a UNIFORM eta grid.
    We compute the differential operator in eta and apply M in x-space.
    Chain form in eta (unknown w(eta)=? ) -- we instead work with s-values at
    x-nodes (eta uniform -> x clustered at 0).  D1_x, D2_x via the chain rule:
       d/dx = (1/x'(eta)) d/deta ;   d2/dx2 = (1/x'^2) d2/deta2 - (x''/x'^3) d/deta
    with x'=p eta^(p-1), x''=p(p-1) eta^(p-2).
    Then A = -(M D1x D2x) with D1x,D2x built on the eta-uniform (x-nonuniform) grid
    but with the CHAIN-RULE coefficients (which are regular in eta!).
    This should avoid the 1/h^2 blowup because the derivative in eta is O(1)."""
    eta = np.linspace(0, 1, N + 2)
    x = eta ** p
    x[0] = 0.0
    h = np.diff(eta)                       # uniform in eta
    # build d/deta and d2/deta2 (standard uniform FD) on interior nodes
    # then D1x_i = (1/xp_i) Deta_i ; D2x_i = (1/xp_i^2) D2eta_i - (xpp_i/xp_i^3) Deta_i
    # assemble full-grid Deta, D2eta for interior nodes i=1..N (eta nodes)
    # Deta (centred, with eta(0)=0 boundary -> s(0)=0 so v_0=0)
    Nint = N
    Deta = np.zeros((Nint, Nint))
    D2eta = np.zeros((Nint, Nint))
    he = h[0]
    for i in range(Nint):
        fi = i + 1
        if i == 0:
            Deta[i, 1] = 0.0
            Deta[i, 0] = 1.0 / (2 * he)         # eta node 1, uses eta0=0 (s=0)
            Deta[i, 1] = 0.0
            # better: centred (v2-v0)/(2h), v0=0
            Deta[i, 0] = 1.0 / (2 * he)
        elif i == Nint - 1:
            Deta[i, i] = 0.0
            Deta[i, i - 1] = -1.0 / (2 * he)
        else:
            Deta[i, i - 1] = -1.0 / (2 * he)
            Deta[i, i + 1] = 1.0 / (2 * he)
        if i == 0:
            D2eta[i, 0] = -2.0 / he ** 2
            D2eta[i, 1] = 1.0 / he ** 2
        elif i == Nint - 1:
            D2eta[i, i - 1] = 2.0 / he ** 2
            D2eta[i, i] = -2.0 / he ** 2
        else:
            D2eta[i, i - 1] = 1.0 / he ** 2
            D2eta[i, i] = -2.0 / he ** 2
            D2eta[i, i + 1] = 1.0 / he ** 2
    etai = eta[1:-1]
    xp = p * etai ** (p - 1)
    xpp = p * (p - 1) * etai ** (p - 2)
    D1x = (1.0 / xp)[:, None] * Deta
    D2x = (1.0 / xp ** 2)[:, None] * D2eta - (xpp / xp ** 3)[:, None] * Deta
    xcols = x[1:-1]
    M = M_phys_fast(xcols, x)
    A = -(M @ D1x @ D2x)
    return xcols, A


if __name__ == "__main__":
    print("chain form in eta with x=eta^p transform")
    for p in [1.0, 1.5, 2.0]:
        print("-- p=%.1f --" % p)
        for N in [40, 80, 160, 320]:
            xc, A = transform_op(N, p)
            w = np.linalg.eigvals(A)
            r = np.sort(-w.real)
            r = r[r > 1e-6]
            print("   N=%3d  slowest rates: %s" % (N, np.round(r[:3], 5)))