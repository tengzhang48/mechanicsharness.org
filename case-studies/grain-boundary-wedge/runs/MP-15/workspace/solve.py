#!/usr/bin/env python3
"""
Stress relaxation by diffusion into a grain boundary of a thin film.

Dimensionless evolution equation (tau = t/t0, x = z/h, sigma normalized by sigma0):
    d sigma/d tau = S[ d^3 sigma/dx^3 ](x),   x in (0,1)
BCs: sigma(0)=0, sigma'(1)=0, sigma''(1)=0;  IC: sigma(x,0)=1.
S is the principal-value integral operator with kernel
    S(x,xi) = 1/(x-xi) - 1/(x+xi) - 2*xi*(x-xi)/(x+xi)^3 .
"""
import json
import sys
from math import comb

import numpy as np
from numpy.polynomial import chebyshev as C
from numpy.polynomial import legendre as L
from scipy.linalg import eig, solve, expm


# ----------------------------------------------------------------------------
# Exact action of S on monomials zeta^m  (used to validate the quadrature)
# ----------------------------------------------------------------------------
def S_monomial_exact(m, z):
    z = float(z)
    P = z**m * np.log(z / (1.0 - z)) - sum(z**(m - 1 - k) / (k + 1) for k in range(m))
    Q = (-z)**m * np.log((1.0 + z) / z) + sum((-z)**(m - 1 - k) / (k + 1) for k in range(m))

    def R(mm):
        out = 0.0
        for j in range(mm + 1):
            c = comb(mm, j) * (-z)**(mm - j)
            Ij = np.log((1.0 + z) / z) if j == 2 else ((1.0 + z)**(j - 2) - z**(j - 2)) / (j - 2)
            out += c * Ij
        return out

    return P - Q - 2.0 * z * R(m + 1) + 2.0 * R(m + 2)


# ----------------------------------------------------------------------------
# Numerical action of S on a list of smooth functions, via split Gauss-Legendre
# ----------------------------------------------------------------------------
def apply_S(zvals, funs, nquad=96, tip_singular=False):
    """Return matrix out[i,k] = S[funs[k]](zvals[i]).  PV handled by subtraction.
    If tip_singular, the true function is funs[k](zeta)/sqrt(1-zeta)."""
    zg, wg = L.leggauss(nquad)
    K = len(funs)
    out = np.zeros((len(zvals), K))
    for i, z in enumerate(zvals):
        b = min(3.0 * z, 1.0)
        intervals = [(0.0, z), (z, b)]
        if b < 1.0 - 1e-15:
            intervals.append((b, 1.0))
        acc = np.zeros(K)
        for (p, q) in intervals:
            if q <= p:
                continue
            if tip_singular and q >= 1.0 - 1e-15:
                smax = np.sqrt(1.0 - p)
                s = 0.5 * smax * (zg + 1.0)
                ze = 1.0 - s * s
                w = 0.5 * smax * wg * 2.0 * s          # dzeta = 2 s ds
            else:
                ze = 0.5 * (p + q) + 0.5 * (q - p) * zg
                w = 0.5 * (q - p) * wg
            F = np.array([f(ze) for f in funs])        # (K, M)
            if tip_singular:
                F = F / np.sqrt(1.0 - ze)
            fz = np.array([f(np.array([z]))[0] for f in funs])
            if tip_singular:
                fz = fz / np.sqrt(1.0 - z)
            I1 = (F - fz[:, None]) / (z - ze)
            I3 = -F / (z + ze)
            I4 = -2.0 * ze * (z - ze) * F / (z + ze)**3
            acc += (I1 + I3 + I4) @ w
        acc += fz * np.log(z / (1.0 - z))
        out[i] = acc
    return out




# ----------------------------------------------------------------------------
# Basis functions with psi(0)=0, d3psi(0)=0, dpsi(1)=0, d2psi(1)=0
# psi_k = T*_{k+3} + a T*_3 + b T*_2 + c T*_1 + d T*_0   (T*_n(z)=T_n(2z-1))
# ----------------------------------------------------------------------------
def _Tn_val(n, p, xend):
    """T_n^{(p)}(xend), xend = +1 or -1 (Chebyshev variable x)."""
    if p > n:
        return 0.0
    if p == 0:
        v = 1.0
    elif p == 1:
        v = float(n * n)
    elif p == 2:
        v = n * n * (n * n - 1) / 3.0
    elif p == 3:
        v = n * n * (n * n - 1) * (n * n - 4) / 15.0
    else:
        raise ValueError
    if xend < 0:
        v *= (-1)**(n + p)
    return v


def build_basis(N):
    basis = []
    ns = [3, 2, 1, 0]
    M = np.array([[ _Tn_val(n, 0, -1) for n in ns],
                  [8.0 * _Tn_val(n, 3, -1) for n in ns],
                  [2.0 * _Tn_val(n, 1,  1) for n in ns],
                  [4.0 * _Tn_val(n, 2,  1) for n in ns]])
    for k in range(1, N + 1):
        n0 = k + 3
        rhs = -np.array([_Tn_val(n0, 0, -1), 8.0 * _Tn_val(n0, 3, -1),
                         2.0 * _Tn_val(n0, 1, 1), 4.0 * _Tn_val(n0, 2, 1)])
        abcd = np.linalg.solve(M, rhs)
        c = np.zeros(n0 + 1)
        c[n0] = 1.0
        for n, v in zip(ns, abcd):
            c[n] += v
        basis.append(c)
    return basis


# ----------------------------------------------------------------------------
# Enriched basis: polynomial part (N funcs) + crack-tip functions (1-x)^{j+1/2},
# j=2..J+1, each corrected by a cubic so all 4 endpoint conditions hold:
#   rho_j(x) = (1-x)^{j+1/2} - 1 - a3*(3x - 3x^2 + x^3),  a3 = -(j+.5)(j-.5)(j-1.5)/6
#   rho_j'''(x) = h_j(x)/sqrt(1-x) - 6*a3,  h_j = -(j+.5)(j-.5)(j-1.5)(1-x)^{j-2}
# ----------------------------------------------------------------------------
def build_basis_enr(N, J):
    items = [('poly', c) for c in build_basis(N)]
    for j in range(2, J + 2):
        items.append(('tip', j))
    return items


def _enr_a3(j):
    return -(j + 0.5) * (j - 0.5) * (j - 1.5) / 6.0


def basis_eval(items, zz):
    zz = np.asarray(zz, float)
    out = np.zeros((len(zz), len(items)))
    for k, (kind, data) in enumerate(items):
        if kind == 'poly':
            out[:, k] = C.chebval(2.0 * zz - 1.0, data)
        else:
            j = data
            a3 = _enr_a3(j)
            out[:, k] = (1.0 - zz)**(j + 0.5) - 1.0 - a3 * (3 * zz - 3 * zz**2 + zz**3)
    return out


def apply_S_basis(z, items, nquad=96):
    """S of (d3/dz3 of each basis function) at points z."""
    smooth_funs, tip_funs = [], []
    tip_smooth_const = []   # the -6*a3 constant accompanying each tip function
    for kind, data in items:
        if kind == 'poly':
            c3 = 8.0 * C.chebder(data, 3)
            smooth_funs.append((True, lambda ze, c3=c3: C.chebval(2.0 * ze - 1.0, c3)))
        else:
            j = data
            a3 = _enr_a3(j)
            fac = -(j + 0.5) * (j - 0.5) * (j - 1.5)
            smooth_funs.append((False, lambda ze, c=-6.0 * a3: c * np.ones_like(ze)))
            tip_funs.append((len(smooth_funs) - 1,
                             lambda ze, j=j, fac=fac: fac * (1.0 - ze)**(j - 2)))
            tip_smooth_const.append(-6.0 * a3)
    # smooth part (tip_singular=False)
    fns = [f for _, f in smooth_funs]
    A = apply_S(z, fns, nquad=nquad, tip_singular=False)
    # tip part (tip_singular=True): add h_j/sqrt(1-x) contributions
    if tip_funs:
        tf = [f for _, f in tip_funs]
        At = apply_S(z, tf, nquad=nquad, tip_singular=True)
        for jj, (col, _) in enumerate(tip_funs):
            A[:, col] += At[:, jj]
    return A


def method1_enr(N, J, nquad=96):
    items = build_basis_enr(N, J)
    M = N + J
    z = np.array([0.5 * (1.0 + np.cos((2 * i - 1) * np.pi / (2 * M)))
                  for i in range(1, M + 1)])
    A = apply_S_basis(z, items, nquad)
    Bm = basis_eval(items, z)
    mu, vecs = eig(A, Bm)
    order = np.argsort(np.abs(mu))
    lam = -mu[order]
    vecs = vecs[:, order]
    v1 = vecs[:, 0]
    f1 = sum(v1[k] * basis_eval([items[k]], np.array([1.0]))[0, 0] for k in range(M))
    fh = sum(v1[k] * basis_eval([items[k]], np.array([0.5]))[0, 0] for k in range(M)) / f1
    return lam, fh, z, items, vecs


def method1(N, nquad=96):
    basis = build_basis(N)
    funs3 = []
    for c in basis:
        c3 = 8.0 * C.chebder(c, 3)          # d3/dz3 = 8 d3/dx3
        funs3.append(lambda ze, c3=c3: C.chebval(2.0 * ze - 1.0, c3))
    z = np.array([0.5 * (1.0 + np.cos((2 * i - 1) * np.pi / (2 * N)))
                  for i in range(1, N + 1)])
    A = apply_S(z, funs3, nquad)
    Bm = np.array([[C.chebval(2.0 * zi - 1.0, c) for c in basis] for zi in z])
    mu, vecs = eig(A, Bm)
    order = np.argsort(np.abs(mu))
    lam = -mu[order]
    vecs = vecs[:, order]
    v1 = vecs[:, 0]
    f1 = sum(v1[k] * np.sum(basis[k]) for k in range(N))   # psi_k(1)=sum of cheb coeffs
    x05 = 2.0 * 0.5 - 1.0
    fh = sum(v1[k] * C.chebval(x05, basis[k]) for k in range(N)) / f1
    return lam, fh, z, basis, vecs


def check_quadrature():
    print("== check: numerical S vs exact monomial integrals ==")
    zs = [0.0137, 0.05, 0.25, 0.5, 0.9, 0.995]
    worst = 0.0
    for m in range(0, 9):
        funs = [lambda ze, m=m: np.asarray(ze)**m]
        num = apply_S(np.array(zs), funs, nquad=128)[:, 0]
        exa = np.array([S_monomial_exact(m, z) for z in zs])
        err = np.max(np.abs(num - exa))
        worst = max(worst, err)
        print(f"  m={m}: max abs err over test z = {err:.3e}")
    print(f"  worst = {worst:.3e}")


# ----------------------------------------------------------------------------
# Method 2: Chebyshev-tau with tail enrichment.
#   basis = {T*_0..T*_{N-1}} + {tail_N((1-x)^{j+1/2}), j=2..J+1}
#   The tail (Chebyshev modes >= N of the bump) has disjoint support from the
#   polynomial block, so the tau matrix stays well conditioned; the span equals
#   polynomials(<N) + bump functions, i.e. captures the crack-tip singularity.
#   Equations: (N+J-4) collocation rows (Gauss pts) + 4 BC rows.
# ----------------------------------------------------------------------------
def _bump(x, j):
    return np.power(1.0 - np.asarray(x, float), j + 0.5)


def _bump_fac(j):
    return -(j + 0.5) * (j - 0.5) * (j - 1.5)   # d3/dx3 (1-x)^{j+1/2} = fac*(1-x)^{j-5/2}


def _bump_tail(j, N, NI=512):
    """Return (tail_cheb, proj_cheb): shifted-Chebyshev coeffs of the bump tail
    (modes >= N) and of its projection (modes < N), via NI-pt interpolation."""
    th = np.pi * (np.arange(NI) + 0.5) / NI
    xg = 0.5 * (1.0 + np.cos(th))
    vg = _bump(xg, j)
    K = np.cos(np.outer(np.arange(NI), th))       # K[k,i] = cos(k th_i)
    a = (2.0 / NI) * (K @ vg)
    a[0] *= 0.5
    tail = a.copy(); tail[:N] = 0.0
    proj = a.copy(); proj[N:] = 0.0
    return tail, proj


def method2_tau(N, J=3, nquad=96, NI=512):
    nbc = 4
    ncol = N + J
    nrow = ncol - nbc
    zc = np.array([0.5 * (1.0 + np.cos((2 * i - 1) * np.pi / (2 * nrow)))
                   for i in range(1, nrow + 1)])
    # ---- B rows: collocation of basis functions
    B = np.zeros((ncol, ncol))
    xx = 2.0 * zc - 1.0
    for j in range(N):
        c = np.zeros(j + 1); c[j] = 1.0
        B[:nrow, j] = C.chebval(xx, c)
    tails, projs, hfuns = [], [], []
    for jj, j in enumerate(range(2, J + 2)):
        tail, proj = _bump_tail(j, N, NI)
        tails.append(tail); projs.append(proj)
        B[:nrow, N + jj] = C.chebval(xx, tail)
        hfuns.append(_bump_fac(j))
    # ---- BC rows in B: phi(0)=0, phi'''(0)=0, phi'(1)=0, phi''(1)=0
    bcs = [(0, 1.0, -1.0), (3, 8.0, -1.0), (1, 2.0, 1.0), (2, 4.0, 1.0)]
    for r, (p, fac, xend) in enumerate(bcs):
        row = np.array([fac * _Tn_val(j, p, xend) for j in range(N)])
        B[nrow + r, :N] = row
        for jj, j in enumerate(range(2, J + 2)):
            # tail = bump - proj ; endpoint functionals:
            bf = _bump_fac(j)
            if r == 0:    # value at 0: bump(0)=1 minus proj(0)
                v = 1.0 - C.chebval(-1.0, projs[jj])
            elif r == 1:  # third deriv at 0: bump'''(0) - proj'''(0)*8
                v = bf - 8.0 * C.chebval(-1.0, C.chebder(projs[jj], 3))
            elif r == 2:  # first deriv at 1: 0 - proj'(1)*2
                v = -2.0 * C.chebval(1.0, C.chebder(projs[jj], 1))
            else:         # second deriv at 1: 0 - proj''(1)*4
                v = -4.0 * C.chebval(1.0, C.chebder(projs[jj], 2))
            B[nrow + r, N + jj] = v
    # ---- A rows: S of third derivatives at collocation pts; zero BC rows
    funs = []
    for j in range(N):
        c = np.zeros(j + 1); c[j] = 1.0
        c3 = 8.0 * C.chebder(c, 3)
        funs.append(lambda ze, c3=c3: C.chebval(2.0 * np.asarray(ze) - 1.0, c3))
    A = np.zeros((ncol, ncol))
    A[:nrow, :N] = apply_S(zc, funs, nquad=nquad, tip_singular=False)
    if J > 0:
        # bump''' = fac*(1-x)^{j-2}/sqrt(1-x)  (tip path) minus projection part
        hf = [lambda ze, f=hfuns[k], j=k + 2: f * np.power(1.0 - np.asarray(ze), j - 2)
              for k in range(J)]
        At = apply_S(zc, hf, nquad=nquad, tip_singular=True)
        pf = []
        for jj in range(J):
            c3 = 8.0 * C.chebder(projs[jj], 3)
            pf.append(lambda ze, c3=c3: C.chebval(2.0 * np.asarray(ze) - 1.0, c3))
        Ap = apply_S(zc, pf, nquad=nquad, tip_singular=False)
        A[:nrow, N:] = At - Ap
    # ---- column scaling of enrichment columns (unit collocation norm)
    for jj in range(J):
        s = np.sqrt(np.sum(B[:nrow, N + jj]**2))
        if s > 0:
            B[:, N + jj] /= s; A[:, N + jj] /= s
    # ---- eliminate BC rows: c = Z y with Z = null(B_bc); reduced pencil
    Bbc = B[nrow:, :].copy()
    for r in range(nbc):                               # row scaling before SVD
        Bbc[r, :] /= np.max(np.abs(Bbc[r, :]))
    _, sv, Vt = np.linalg.svd(Bbc)
    Z = Vt[nbc:].T                                     # (ncol) x (ncol-4)
    Ar = A[:nrow, :] @ Z
    Br = B[:nrow, :] @ Z
    mu, vecs = eig(Ar, Br)
    order = np.argsort(np.abs(mu))
    lam = -mu[order]
    vecs = Z @ vecs[:, order]                          # back to basis coeffs
    # normalize mode 1 to f(1)=1 and evaluate f(1/2)
    def modeval(v, x):
        xx = 2.0 * x - 1.0
        val = 0.0
        for j in range(N):
            c = np.zeros(j + 1); c[j] = 1.0
            val = val + v[j] * C.chebval(xx, c)
        for jj in range(J):
            val = val + v[N + jj] * C.chebval(xx, tails[jj])
        return val
    v1 = vecs[:, 0]
    fh = modeval(v1, 0.5) / modeval(v1, 1.0)
    return lam, fh, (B, A, tails)


# ----------------------------------------------------------------------------
# Method 1t: spectral collocation in the mapped variable t, x = 2t - t^2
# (s = sqrt(1-x) = 1-t).  In s the crack-tip series is polynomial; the four
# physical endpoint conditions at the tip become  c1=c2=c3=c4=0 in the
# s-power-series, i.e. psi'(1)=psi''(1)=psi'''(1)=psi''''(1)=0; at the mouth
# psi(0)=0 and psi'''(0)+psi''(0)+psi'(0)=0 (zero x-third-derivative at x=0).
# For phi(t) = sum c_n (1-t)^n,  d3phi/dx3 = H(s)/s with
# H(s) = -sum_{n>=5} c_n n(n-2)(n-4)/8 * s^(n-5)  (exact, well conditioned).
# ----------------------------------------------------------------------------
def _mono_fun(n, p, tend):
    """(d/dt)^p of t^n at t=tend."""
    if n < p:
        return 0.0
    v = 1.0
    for i in range(p):
        v *= (n - i)
    return v * tend**(n - p)


def _Tstar_fun(m, p, tend):
    return (2.0**p) * _Tn_val(m, p, 2.0 * tend - 1.0)


def _fun6(p):
    """The 6 endpoint functionals of a t-polynomial p (power-basis coeffs)."""
    F1 = np.polyval(p, 0.0)
    p1 = np.polyder(p, 1); p2 = np.polyder(p, 2); p3 = np.polyder(p, 3)
    p4 = np.polyder(p, 4)
    F2 = np.polyval(p3, 0.0) + np.polyval(p2, 0.0) + np.polyval(p1, 0.0)
    F3 = np.polyval(p1, 1.0); F4 = np.polyval(p2, 1.0); F5 = np.polyval(p3, 1.0)
    F6 = np.polyval(p4, 1.0)
    return np.array([F1, F2, F3, F4, F5, F6])


def build_basis_t(N):
    """psi_k(t) = T*_{k+5}(t) + cubic-... degree<=5 correction (power basis)."""
    M6 = np.array([[_mono_fun(n, 0, 0.0),
                    _mono_fun(n, 3, 0.0) + _mono_fun(n, 2, 0.0) + _mono_fun(n, 1, 0.0),
                    _mono_fun(n, 1, 1.0), _mono_fun(n, 2, 1.0),
                    _mono_fun(n, 3, 1.0), _mono_fun(n, 4, 1.0)] for n in range(6)]).T
    basis = []
    for k in range(1, N + 1):
        m = k + 5
        rhs = -np.array([_Tstar_fun(m, 0, 0.0),
                         _Tstar_fun(m, 3, 0.0) + _Tstar_fun(m, 2, 0.0) + _Tstar_fun(m, 1, 0.0),
                         _Tstar_fun(m, 1, 1.0), _Tstar_fun(m, 2, 1.0),
                         _Tstar_fun(m, 3, 1.0), _Tstar_fun(m, 4, 1.0)])
        d = np.linalg.solve(M6, rhs)
        # power-basis coeffs of T*_m (ascending order)
        cc = np.zeros(m + 1)
        cc[m] = 1.0
        tc = C.cheb2poly(cc)          # coeffs of T_m(x), ascending in x
        pw = np.array([tc[j] * (2.0**j) * math.comb(0, 0) for j in range(m + 1)])
        # T*_m(t) = T_m(2t-1): expand binomially
        pw = np.zeros(m + 1)
        for j in range(m + 1):
            for l in range(j + 1):
                pw[l] += tc[j] * math.comb(j, l) * (2.0**l) * ((-1.0)**(j - l))
        pw = C.poly2cheb(pw)  # placeholder no-op guard (unused)
        basis.append(None)
    return basis, M6


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("quad", "all"):
        check_quadrature()
    if stage in ("eig", "all"):
        print("== method 1 (spectral collocation) convergence ==")
        for N in [12, 16, 20, 24, 28, 32, 40]:
            lam, fh, _, _, _ = method1(N)
            print(f" N={N:3d}  lambda1={lam[0].real:.10f} (imag {lam[0].imag:.2e})"
                  f"  lambda2={lam[1].real:.10f} (imag {lam[1].imag:.2e})"
                  f"  f(0.5)={fh.real:.10f} (imag {fh.imag:.2e})")

