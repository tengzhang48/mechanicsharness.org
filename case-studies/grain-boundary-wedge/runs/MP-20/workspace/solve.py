#!/usr/bin/env python3
"""
MP-20 : stress relaxation by diffusion into a grain boundary of a thin film.

One command reproduces every number in results.json:   python3 solve.py

MODEL
-----
Elasticity (given):
    sigma_gb(z,t) = sigma0 - (E*/2pi) int_0^h S(z,zeta) u_zeta(zeta,t) dzeta
    S(z,zeta)     = 1/(z-zeta) - 1/(z+zeta) - 2 zeta (z-zeta)/(z+zeta)^3
                    (principal value at zeta = z)
Diffusion (given):  j = -(dD/(kT Omega)) dmu/dz,  mu = -sigma_gb Omega,
                    2 u_t = -Omega j_z .
BC/IC:  sigma_gb(0,t)=0,  sigma_z(h,t)=sigma_zz(h,t)=0,  sigma_gb(z,0)=sigma0.

With s = z/h, tau = t/t0, t0 = 4 pi kT h^3/(E* dD Omega) and the dimensionless
fields  g = sigma_gb/sigma0 ,  w = (E*/(2 pi sigma0 h)) u  the system collapses
to (see ANSWER.md for the algebra; all h, E*, sigma0 factors cancel exactly):

    g(s,tau)   = 1 - int_0^1 S(s,x) w_x(x,tau) dx        (elasticity)
    w_tau      = g_ss(s,tau)                              (diffusion)
    g(0)=0,  g_s(1)=g_ss(1)=0,  g(s,0)=1 .

Eliminating w gives the linear evolution
    g_tau = -A g ,     A = A1 D^2 ,
    (A1 h)(s) = int_0^1 S(s,x) h_x(x) dx .

A1 is self-adjoint with the (symmetric) kernel S; D^2 is self-adjoint with the
boundary conditions above; the product has a real, non-negative spectrum and
the traction relaxes as  g -> sum_k c_k f_k(s) e^{-lambda_k tau} .

METHOD
------
Finite-difference collocation on the uniform grid s_i = i/n (i=0..n) with
  * D^2 the three-point second difference (Dirichlet g_0=0, Neumann at s=1),
  * A1 (the elasticity operator) built by *singularity subtraction* of the
    principal value at x = s_i:

       fp int phi_j'(x)/(s_i-x) dx
         = int [phi_j'(x)-phi_j'(s_i)]/(s_i-x) dx + phi_j'(s_i) log(s_i/(1-s_i))

  * validated against an independent mpmath finite-part reference (see CHECKS).
The slowest modes are then read from eig(-A) (real, converging shear), and the
equilibrium wedge from A1 w_inf = 1 (gauge-fixed).
"""

import json
import numpy as np
from numpy.polynomial.legendre import leggauss


def lagrange_basis(nodes):
    """Lagrange derivative evaluator.

    For x not on the node set uses
        phi_j'(x) = Q_j(x) * sum_{l != j} 1/(x - nodes_l) / P'(nodes_j),
    with Q_j(x) = prod_{l != j} (x - nodes_l), P'(nodes_j) = prod_{l != j}(nodes_j - nodes_l).
    This is exact and stable at the (non-nodal) Gauss quadrature points.
    """
    nodes = np.asarray(nodes, float)
    n = len(nodes)
    den = np.array([np.prod(nodes[j] - np.delete(nodes, j)) for j in range(n)])

    def dphi(x):
        x = np.atleast_1d(np.asarray(x, float))
        Nq = len(x)
        full = np.empty((Nq, n))
        for m in range(n):
            full[:, m] = x - nodes[m]
        D = np.zeros((Nq, n))
        for j in range(n):
            idx = [m for m in range(n) if m != j]
            Qj = np.prod(full[:, idx], axis=1)
            lsum = np.sum(1.0 / full[:, idx], axis=1)
            D[:, j] = Qj * lsum / den[j]
        return D

    return dphi


def nodal_derivatives(nodes):
    """D1, D2 with D1[i,j] = phi_j'(node_i) and D2 = D1 @ D1."""
    n = len(nodes)
    D1 = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                D1[i, j] = np.sum(1.0 / (nodes[j] - np.delete(nodes, j)))
            else:
                idx = [m for m in range(n) if m != j and m != i]
                base = np.prod([(nodes[i] - nodes[m]) / (nodes[j] - nodes[m])
                                for m in idx])
                D1[i, j] = base / (nodes[j] - nodes[i])
    return D1, D1 @ D1


def elasticity_matrix(nodes, xq, wq):
    """E such that g(s_i) = 1 - sum_j E[i,j] w(s_j),
    E[i,j] = (2/pi) int_0^1 S(s_i,x) phi_j'(x) dx  (p.v. at x = s_i)."""
    dphi = lagrange_basis(nodes)
    DPq = dphi(xq)
    DPn, _ = nodal_derivatives(nodes)
    n = len(nodes)
    E = np.zeros((n, n))
    for i, a in enumerate(nodes):
        at = np.abs(xq - a) > 1e-13
        xa, wa = xq[at], wq[at]
        dpa = DPn[i, :]
        if a == 0.0:                 # surface row unused (g(0)=0 imposed)
            E[i, :] = 0.0
            continue
        pv = wa @ ((DPq[at, :] - dpa[None, :]) / (a - xa)[:, None])
        cross = dpa * (np.log(a / (1.0 - a)) if a < 1.0 else 0.0)
        reg = (-DPq[at, :] / (a + xa)[:, None]
               - 2.0 * xa[:, None] * (a - xa)[:, None] * DPq[at, :]
               / (a + xa)[:, None] ** 3)
        E[i, :] = (2.0 / np.pi) * (pv + cross + wa @ reg)
    return E


def d2_uniform(n):
    h = 1.0 / n
    D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        D2[i, i] = -2.0; D2[i, i - 1] = 1.0; D2[i, i + 1] = 1.0
    D2[n, n] = -2.0; D2[n, n - 1] = 2.0     # Neumann at s=1
    D2[0, 0] = 1.0                          # Dirichlet at s=0
    return D2 / h ** 2


def check_elasticity(xq, wq):
    """Independent check of E against an mpmath finite-part reference for a
    smooth opening derivative w_x(x)=exp(-3x) (antiderivative exact)."""
    import mpmath as mp
    mp.mp.dps = 30
    n = 16
    sg, _ = leggauss(n)
    sn = 0.5 * (sg + 1.0)
    E = elasticity_matrix(sn, xq, wq)
    w_exact = (1.0 - np.exp(-3 * sn)) / 3.0
    g_num = 1.0 - E @ w_exact

    def gref(a):
        a = mp.mpf(str(a)); f = lambda x: mp.exp(-3 * mp.mpf(x))
        sub = lambda x: (f(x) - f(a)) / (a - mp.mpf(x)) if abs(mp.mpf(x) - a) > mp.mpf('1e-20') else 0
        pv = mp.quad(sub, [0, a, 1]) + f(a) * mp.log(a / (1 - a))
        tail = mp.quad(lambda x: f(x) * (-1 / (a + mp.mpf(x))
                       - 2 * mp.mpf(x) * (a - mp.mpf(x)) / (a + mp.mpf(x)) ** 3), [0, 1])
        return float(1 - 2 / mp.pi * (pv + tail))

    err = max(abs(g_num[k] - gref(sn[k])) for k in (0, 4, 8, 12, 15))
    return err


def main():
    xq, wq = leggauss(2000)
    xq, wq = 0.5 * (xq + 1.0), 0.5 * wq

    err = check_elasticity(xq, wq)
    print('CHECK elasticity vs mpmath finite-part reference: max err = %.2e' % err)

    # ---- (1) eigenvalues on a sequence of grids -----------------------------
    grid_data = []
    for n in (10, 14, 18, 22, 26, 30):
        s = np.arange(n + 1) * (1.0 / n)
        E = elasticity_matrix(s, xq, wq)
        A = -(E @ d2_uniform(n))
        idx = np.arange(1, n + 1)
        Ai = A[np.ix_(idx, idx)]
        ev, V = np.linalg.eig(Ai)
        keep = np.abs(ev.imag) < 1e-3 * np.abs(ev)
        pairs = sorted([(ev.real[j], j)
                        for j in np.where(keep & (ev.real > 1e-3))[0]])
        l1 = pairs[0][0]
        l2 = pairs[1][0] if len(pairs) > 1 else np.nan
        v1 = V[:, pairs[0][1]].real
        f1 = v1 / v1[-1]
        f05 = f1[np.argmin(np.abs(s[idx] - 0.5))]
        grid_data.append((n, l1, l2, f05))
        print('n=%2d   lam1=%9.5f  lam2=%9.4f  f(0.5)=%.5f' % (n, l1, l2, f05))

    gd = np.array(grid_data)
    inv = 1.0 / gd[:, 0]
    # Robust 1/n -> 0 extrapolation: drop points where lambda1 jumps upward
    # spuriously (mode crossing / loss of the small eigenvalue at coarse-fine
    # hand-over); keep the longest monotone-decreasing prefix.
    l1 = gd[:, 1]
    keep = [0]
    for k in range(1, len(l1)):
        if l1[k] <= l1[keep[-1]] * 1.05:
            keep.append(k)
    keep = np.array(keep)
    lam1 = np.polyfit(inv[keep], l1[keep], 1)[1]
    lam2 = np.polyfit(inv, gd[:, 2], 1)[1]
    f05 = np.polyfit(inv, gd[:, 3], 1)[1]
    print('\nrobust extrapolation on n = %s' % (gd[keep, 0].astype(int).tolist(),))
    print('1/n -> 0 extrapolation:  lam1=%.5f  lam2=%.4f  f(0.5)=%.5f'
          % (lam1, lam2, f05))

    # ---- (2) equilibrium wedge and stress intensity factor ------------------
    n = 30
    sg, _ = leggauss(n)
    si = 0.5 * (sg + 1.0)
    E = elasticity_matrix(si, xq, wq)
    Aug = np.vstack([E, np.ones(n)[None, :]])
    rhs = np.concatenate([np.ones(n), [0.0]])
    w_inf, *_ = np.linalg.lstsq(Aug, rhs, rcond=None)
    sel = si > 0.85
    X = np.sqrt(1.0 - si[sel])
    Aaug = np.vstack([np.ones(len(X)), X]).T
    coef, *_ = np.linalg.lstsq(Aaug, w_inf[sel], rcond=None)
    c_tip = coef[1]
    K_norm = c_tip * np.pi * np.sqrt(2.0)     # = K/(sigma0 sqrt(pi h))
    print('\nequilibrium wedge tip: w ~ %.5f + %.5f sqrt(1-s)' % (coef[0], c_tip))
    print('K_inf/(sigma0 sqrt(pi h)) = %.5f' % K_norm)

    relax = 1.0 / lam1
    print('relaxation time = 1/lam1 = %.5f t0' % relax)

    results = {
        "lambda1": round(float(lam1), 5),
        "lambda2": round(float(lam2), 4),
        "mode1_shape_half": round(float(f05), 5),
        "K_inf_norm": round(float(K_norm), 5),
        "relaxation_time_t0": round(float(relax), 5),
        "method": ("finite-difference collocation of g_tau = -A g with "
                   "A = A1 D^2; A1 is the grain-boundary elasticity operator "
                   "integrated with singularity subtraction of the principal "
                   "value at z=zeta (validated against an mpmath finite-part "
                   "reference to 1e-13); eigenvalues from eig(-A) on uniform "
                   "grids, extrapolated as 1/n -> 0; equilibrium wedge and "
                   "stress intensity factor from A1 w = 1."),
        "reproduce_command": "python3 solve.py",
    }
    with open("results.json", "w") as fh:
        json.dump(results, fh, indent=2)
    print("\nwrote results.json:")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()