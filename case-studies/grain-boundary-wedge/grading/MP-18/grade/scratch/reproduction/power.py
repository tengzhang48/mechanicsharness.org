import numpy as np
import scipy.linalg as sla
from fd_nonuniform import build_D1D2
from phys_oper import M_phys


def march_s(nodes, T=30.0, nsteps=3000):
    """March the s-equation.  Unknowns: s at interior nodes (s(0)=0).
    Reduced evolution  s_tau = -M D1 D2 s  is a pure decay; we instead march the
    physical chain  w_tau = D1 D2 s ,  s = 1 - M w  with the IC s=const, using the
    consistent operator A = -M D1 D2 on the *deviation* d = s - s_ss where s_ss
    is the steady state (s_tau=0 => D1 D2 s_ss = 0).  For a decay test we simply
    iterate  s_{new} = (I - dt A) s_old  implicitly and watch ||s||.
    """
    xcols = np.asarray(nodes[1:-1], float)
    n = len(xcols)
    D1, D2 = build_D1D2(nodes)
    M = M_phys(xcols, nodes)
    A = -(M @ D1 @ D2)               # s_tau = A s
    dt = T / nsteps
    # backward Euler on the homogeneous decay of the deviation:
    # we want the *slowest* mode; use the power method on (I - dt A)^{-1} (which
    # has spectral radius 1 - dt*lambda_smallest -> 1 for the slowest mode).
    Bm = np.eye(n) - dt * A
    lu = sla.lu_factor(Bm)
    s = np.ones(n)
    s = s / np.linalg.norm(s)
    ratios = []
    for k in range(nsteps):
        s_new = sla.lu_solve(lu, s)
        nrm = np.linalg.norm(s_new)
        s_new = s_new / nrm
        # Rayleigh: growth factor per step ~ 1/(1 - dt*(-lambda)) ...
        ratios.append(nrm)
        s = s_new
    # the dominant eigenvalue of Bm^{-1} is 1/(1 - dt*mu), mu largest A eigenvalue
    dom = ratios[-1]
    mu = (1.0 / dom - 1.0) / (-dt) if dom != 0 else np.nan
    return s, mu, ratios


if __name__ == "__main__":
    for N in [40, 80]:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        s, mu, ratios = march_s(nodes, T=30.0, nsteps=3000)
        print("N=%3d  power-method dominant rate (should be slowest) = %.6f" % (N, -mu))
        v = s / s[-1]
        print("   shape f(0.1,0.5,1)=%s" % np.round(np.interp([0.1, 0.5, 1.0], x, v), 4))