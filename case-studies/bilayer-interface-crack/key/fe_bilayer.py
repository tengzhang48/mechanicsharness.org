"""Finite-element solution of a semi-infinite interface crack between two elastic layers (plane strain), written for the
answer key of the bilayer benchmark. Independent of the source paper's integral-equation method.

Geometry (lengths in units of the upper-layer thickness h = 1): layer 1 occupies 0 < x2 < h, layer 2 occupies
-H < x2 < 0. The interface x2 = 0 is cracked for x1 < 0 (traction-free faces) and bonded for x1 > 0. The model is
truncated at x1 = -Lc and x1 = +Lb (both many thicknesses from the tip).

Loading: far behind the tip the arms carry beam stresses (exact for a homogeneous strip under end force and moment):
  upper arm:  sigma11 = N/h - 12 M (x2 - h/2)/h^3            (M > 0 opens the crack)
  lower arm:  sigma11 = -N/H + 12 (M - N (h + H)/2) (x2 + H/2)/H^3
which makes the bonded region far ahead stress-free. They are applied as tractions on x1 = -Lc; x1 = +Lb is free.

Elements: 9-node Lagrange quadrilaterals on a tensor-product grid graded geometrically towards the tip and the interface;
3 x 3 Gauss points. Crack-face nodes are duplicated for x1 < 0. Rigid motion removed by three constraints at x1 = Lb.

Outputs (function solve): the complex interface stress intensity factor K = K1 + i K2 defined by
  sigma22 + i sigma12 = K r^(i eps)/sqrt(2 pi r) on the interface ahead of the tip,
extracted from the crack opening (delta2 + i delta1 = (c1 + c2) K sqrt(r) r^(i eps)/(2 sqrt(2 pi) (1 + 2 i eps) cosh(pi eps)),
delta = upper face minus lower face), extrapolated to r -> 0; and J from a domain integral (G = J).
Run as a script for a quick self-check: python3 -I fe_bilayer.py
"""
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def graded(d0, g, dmax, length, far=0.0):
    """Points 0 = p0 < p1 < ... < pn = length: steps d0, d0 g, ... capped at max(dmax, far * distance), last step
    stretched to fit (far = 0: a fixed cap)."""
    pts, s, d = [0.0], 0.0, d0
    while s + d < length - 1e-12:
        s += d
        pts.append(s)
        d = min(d * g, max(dmax, far * s))
    if length - pts[-1] < 0.3 * d and len(pts) > 1:
        pts[-1] = length
    else:
        pts.append(length)
    return np.array(pts)


def plane_strain_D(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    return np.array([[lam + 2 * mu, lam, 0], [lam, lam + 2 * mu, 0], [0, 0, mu]])


def dundurs(E1, nu1, E2, nu2):
    mu1, mu2 = E1 / (2 * (1 + nu1)), E2 / (2 * (1 + nu2))
    k1, k2 = 3 - 4 * nu1, 3 - 4 * nu2
    G = mu1 / mu2
    alpha = (G * (k2 + 1) - (k1 + 1)) / (G * (k2 + 1) + (k1 + 1))
    beta = (G * (k2 - 1) - (k1 - 1)) / (G * (k2 + 1) + (k1 + 1))
    eps = math.log((1 - beta) / (1 + beta)) / (2 * math.pi)
    c1, c2 = (k1 + 1) / mu1, (k2 + 1) / mu2
    return alpha, beta, eps, c1, c2


# 1D quadratic shape functions on [-1, 1] at nodes -1, 0, 1
def n1d(t):
    return np.array([0.5 * t * (t - 1), 1 - t * t, 0.5 * t * (t + 1)])


def dn1d(t):
    return np.array([t - 0.5, -2 * t, t + 0.5])


GP = np.array([-math.sqrt(0.6), 0.0, math.sqrt(0.6)])
GW = np.array([5 / 9, 8 / 9, 5 / 9])


def beam_G(E1, nu1, E2, nu2, h, H, N, M):
    """Energy release rate from the energy far behind the tip (beam stresses; far ahead is stress-free)."""
    Eb1, Eb2 = E1 / (1 - nu1 ** 2), E2 / (1 - nu2 ** 2)
    M2 = M - N * (h + H) / 2.0          # moment of the lower arm about its mid-plane: c H^3/12
    return N ** 2 / (2 * Eb1 * h) + 6 * M ** 2 / (Eb1 * h ** 3) + N ** 2 / (2 * Eb2 * H) + 6 * M2 ** 2 / (Eb2 * H ** 3)


def solve(E1, nu1, E2, nu2, H, N, M, h=1.0, d0=2e-4, g=1.15, Lfac=8.0, fit=(4e-3, 4e-2), refine=1.0, far=0.0, verbose=False):
    T = h + H
    Lc = Lb = Lfac * T
    d0 = d0 / refine
    g = g ** (1.0 / refine)
    # grid lines (element corners)
    xr = graded(d0, g, 0.25 * min(h, H) / refine, Lb, far / refine)
    xl = graded(d0, g, 0.25 * min(h, H) / refine, Lc, far / refine)
    xc = np.concatenate([-xl[::-1], xr[1:]])
    yu = graded(d0, g, 0.06 * h / refine, h)
    yd = graded(d0, g, (0.06 * H if far == 0 else 0.06 * h) / refine, H, far / refine)
    yc = np.concatenate([-yd[::-1], yu[1:]])
    # node grid: add midpoints
    def with_mid(c):
        out = np.empty(2 * len(c) - 1)
        out[0::2] = c
        out[1::2] = 0.5 * (c[:-1] + c[1:])
        return out
    X = with_mid(xc)
    Y = with_mid(yc)
    nx, ny = len(X), len(Y)
    j0 = 2 * (len(yd) - 1)                      # node row of the interface
    i0 = 2 * (len(xl) - 1)                      # node column of the tip
    assert abs(Y[j0]) < 1e-14 and abs(X[i0]) < 1e-14
    nid = np.arange(nx * ny).reshape(ny, nx)    # nid[j, i]
    # duplicate nodes on the crack faces (x1 < 0) for the lower layer
    dup = {}
    nn = nx * ny
    for i in range(i0):
        dup[i] = nn
        nn += 1
    coords = np.zeros((nn, 2))
    coords[:nx * ny, 0] = np.tile(X, ny)
    coords[:nx * ny, 1] = np.repeat(Y, nx)
    for i, k in dup.items():
        coords[k] = (X[i], 0.0)
    # elements
    ne_x, ne_y = len(xc) - 1, len(yc) - 1
    conn = np.zeros((ne_x * ne_y, 9), dtype=np.int64)
    mat = np.zeros(ne_x * ne_y, dtype=np.int64)
    e = 0
    for ey in range(ne_y):
        for ex in range(ne_x):
            jj, ii = 2 * ey, 2 * ex
            nodes = []
            for b in range(3):          # local ordering: (a, b) -> a + 3 b, a along x1, b along x2
                for a in range(3):
                    j, i = jj + b, ii + a
                    k = nid[j, i]
                    if j == j0 and i < i0 and ey < len(yd) - 1:   # element below the interface touching the crack
                        k = dup[i]
                    nodes.append(k)
            conn[e] = nodes
            mat[e] = 1 if ey >= len(yd) - 1 else 2
            e += 1
    # element sizes
    ex_idx = np.tile(np.arange(ne_x), ne_y)
    ey_idx = np.repeat(np.arange(ne_y), ne_x)
    hx = (xc[1:] - xc[:-1])[ex_idx]
    hy = (yc[1:] - yc[:-1])[ey_idx]
    D = {1: plane_strain_D(E1, nu1), 2: plane_strain_D(E2, nu2)}
    # stiffness, vectorized over elements: B for each Gauss point (rectangles: dx = hx/2 dt)
    rows, cols, vals = [], [], []
    Ke = np.zeros((len(conn), 18, 18))
    for p, (tp, wp) in enumerate(zip(GP, GW)):
        for q, (tq, wq) in enumerate(zip(GP, GW)):
            Na, Nb = n1d(tp), n1d(tq)
            dNa, dNb = dn1d(tp), dn1d(tq)
            dNdx_ref = np.array([dNa[a] * Nb[b] for b in range(3) for a in range(3)])   # d/dt
            dNdy_ref = np.array([Na[a] * dNb[b] for b in range(3) for a in range(3)])   # d/du
            dNdx = dNdx_ref[None, :] * (2 / hx)[:, None]
            dNdy = dNdy_ref[None, :] * (2 / hy)[:, None]
            B = np.zeros((len(conn), 3, 18))
            B[:, 0, 0::2] = dNdx
            B[:, 1, 1::2] = dNdy
            B[:, 2, 0::2] = dNdy
            B[:, 2, 1::2] = dNdx
            detJ = hx * hy / 4
            for m in (1, 2):
                sel = mat == m
                Ke[sel] += np.einsum('eki,kl,elj->eij', B[sel], D[m], B[sel]) * (wp * wq * detJ[sel])[:, None, None]
    dofs = np.zeros((len(conn), 18), dtype=np.int64)
    dofs[:, 0::2] = 2 * conn
    dofs[:, 1::2] = 2 * conn + 1
    I = np.repeat(dofs, 18, axis=1).ravel()
    Jx = np.tile(dofs, (1, 18)).ravel()
    K = sp.csr_matrix((Ke.ravel(), (I, Jx)), shape=(2 * nn, 2 * nn))
    # loads on x1 = -Lc: traction t = sigma . n with n = -e1 -> t1 = -sigma11
    F = np.zeros(2 * nn)
    def s_upper(y):
        return N / h - 12 * M * (y - h / 2) / h ** 3
    def s_lower(y):
        return -N / H + 12 * (M - N * (h + H) / 2) * (y + H / 2) / H ** 3
    for ey in range(ne_y):
        jj = 2 * ey
        nodes = []
        for b in range(3):
            j = jj + b
            k = nid[j, 0]
            if j == j0 and ey < len(yd) - 1:
                k = dup[0]
            nodes.append(k)
        y0, y1 = yc[ey], yc[ey + 1]
        sfun = s_upper if ey >= len(yd) - 1 else s_lower
        for tq, wq in zip(GP, GW):
            y = 0.5 * (y0 + y1) + 0.5 * (y1 - y0) * tq
            Nb = n1d(tq)
            for b in range(3):
                F[2 * nodes[b]] += -sfun(y) * Nb[b] * wq * 0.5 * (y1 - y0)
    # constraints at x1 = Lb: node (Lb, 0) fixed, node (Lb, h) u1 = 0
    fixed = [2 * nid[j0, nx - 1], 2 * nid[j0, nx - 1] + 1, 2 * nid[ny - 1, nx - 1]]
    free = np.setdiff1d(np.arange(2 * nn), fixed)
    u = np.zeros(2 * nn)
    Kff = K[free][:, free].tocsc()
    u[free] = spla.spsolve(Kff, F[free], permc_spec='COLAMD')
    # crack opening on the faces
    alpha, beta, eps, c1, c2 = dundurs(E1, nu1, E2, nu2)
    ii = np.arange(i0 - 1, -1, -1)                       # nodes behind the tip, nearest first
    r = -X[ii]
    up = np.array([nid[j0, i] for i in ii])
    lo = np.array([dup[i] for i in ii])
    d1 = u[2 * up] - u[2 * lo]
    d2 = u[2 * up + 1] - u[2 * lo + 1]
    fac = (c1 + c2) / (2 * math.sqrt(2 * math.pi) * (1 + 2j * eps) * math.cosh(math.pi * eps))
    Kapp = (d2 + 1j * d1) / (fac * np.sqrt(r) * r ** (1j * eps))
    sel = (r >= fit[0] * h) & (r <= fit[1] * h)
    A = np.vstack([np.ones(sel.sum()), r[sel], r[sel] ** 2]).T
    coef = np.linalg.lstsq(A, Kapp[sel], rcond=None)[0]
    Kc = coef[0]
    # J by the domain integral, q = 1 for rho < R1, linear to 0 at R2
    R1, R2 = 0.1 * min(h, H), 0.3 * min(h, H)
    J = 0.0
    for p, (tp, wp) in enumerate(zip(GP, GW)):
        for q_, (tq, wq) in enumerate(zip(GP, GW)):
            Na, Nb = n1d(tp), n1d(tq)
            dNa, dNb = dn1d(tp), dn1d(tq)
            Nv = np.array([Na[a] * Nb[b] for b in range(3) for a in range(3)])
            dNdx = np.array([dNa[a] * Nb[b] for b in range(3) for a in range(3)])[None, :] * (2 / hx)[:, None]
            dNdy = np.array([Na[a] * dNb[b] for b in range(3) for a in range(3)])[None, :] * (2 / hy)[:, None]
            xe = coords[conn, 0]
            ye = coords[conn, 1]
            xg = xe @ Nv
            yg = ye @ Nv
            rho = np.hypot(xg, yg)
            near = rho < R2 * 1.5
            if not near.any():
                continue
            ue1 = u[2 * conn]
            ue2 = u[2 * conn + 1]
            u11 = np.sum(dNdx * ue1, 1); u12 = np.sum(dNdy * ue1, 1)
            u21 = np.sum(dNdx * ue2, 1); u22 = np.sum(dNdy * ue2, 1)
            # q gradient from nodal q
            qn = np.clip((R2 - np.hypot(xe, ye)) / (R2 - R1), 0, 1)
            q1 = np.sum(dNdx * qn, 1); q2 = np.sum(dNdy * qn, 1)
            for m in (1, 2):
                Dm = D[m]
                s11 = Dm[0, 0] * u11 + Dm[0, 1] * u22
                s22 = Dm[1, 0] * u11 + Dm[1, 1] * u22
                s12 = Dm[2, 2] * (u12 + u21)
                W = 0.5 * (s11 * u11 + s22 * u22 + s12 * (u12 + u21))
                integrand = (s11 * u11 + s12 * u21 - W) * q1 + (s12 * u11 + s22 * u21) * q2
                sel_m = (mat == m) & near
                J += np.sum(integrand[sel_m] * (hx * hy / 4)[sel_m]) * wp * wq
    Kabs_J = math.sqrt(16 * math.cosh(math.pi * eps) ** 2 * J / (c1 + c2))
    out = {'alpha': alpha, 'beta': beta, 'eps': eps, 'K': Kc, 'J': J, 'G_beam': beam_G(E1, nu1, E2, nu2, h, H, N, M),
           'Kabs_J': Kabs_J, 'psi_deg': math.degrees(np.angle(Kc * h ** (1j * eps))), 'ndof': 2 * nn,
           'r': r, 'Kapp': Kapp}
    if verbose:
        print('ndof %d  alpha %.4f beta %.4f eps %.5f  J %.6g  G_beam %.6g  |K| cod %.6g  |K| J %.6g  psi %.4f deg' %
              (2 * nn, alpha, beta, eps, J, out['G_beam'], abs(Kc), Kabs_J, out['psi_deg']))
    return out


if __name__ == '__main__':
    import time
    t0 = time.time()
    # homogeneous, equal thicknesses: exact K_I = sqrt(3) P h^-1/2 + 2 sqrt(3) M h^-3/2, K_II = 2 P h^-1/2 (paper Eq. 2.19),
    # with P the compressive force in the upper arm. solve() takes the tensile force N = -P.
    for (Nn, Mm) in ((1.0, 0.0), (0.0, 1.0)):
        o = solve(1.0, 0.3, 1.0, 0.3, 1.0, Nn, Mm, verbose=True)
        print('   exact K = %s' % (complex(-math.sqrt(3) * Nn + 2 * math.sqrt(3) * Mm, -2 * Nn)), ' FE K = %.5f%+.5fi' % (o['K'].real, o['K'].imag),
              ' time %.0f s' % (time.time() - t0))
