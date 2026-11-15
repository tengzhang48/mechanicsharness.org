"""Check of the bimaterial part of fe_bilayer.py against an exact solution that does not involve layers: a crack of
length 2a on the interface of two bonded half-planes under remote sigma22 = sigma, sigma12 = tau. Right tip:
  K = (sigma + i tau) (1 + 2 i eps) sqrt(pi a) (2a)^(-i eps)        (J. R. Rice 1988, J. Appl. Mech. 55, 98-103),
with sigma22 + i sigma12 = K r^(i eps)/sqrt(2 pi r) ahead of the tip, material 1 above. This tests the oscillatory
singularity, the sign of eps, the crack-opening formula and the extraction used by fe_bilayer.py.
The half-planes are a square of side 2W = 80a with the uniform remote stress applied as tractions (sigma11 in material 1
chosen so that the remote strain eps11 is continuous). Both tips are refined. Element routines from fe_bilayer.py.
Run: python3 -I check_griffith_interface.py   (NumPy, SciPy). Exit 1 if the phase differs by more than 0.1 degree.
"""
import math, os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fe_bilayer import graded, plane_strain_D, dundurs, n1d, dn1d, GP, GW  # noqa: E402


def run(E1, nu1, E2, nu2, sig, tau, a=1.0, W=40.0, d0=2e-4, g=1.15):
    # x1 grid: crack on (-2a, 0); refine at both tips
    out_r = graded(d0, g, 2.0, W - 0.0)                     # from the right tip outwards
    inner = graded(d0, g, 0.1, a)                           # from each tip to the crack centre
    out_l = graded(d0, g, 2.0, W - 2 * a)
    xc = np.unique(np.round(np.concatenate([-2 * a - out_l[::-1], -2 * a + inner, -inner[::-1], out_r]), 14))
    ys = graded(d0, g, 2.0, W)
    yc = np.concatenate([-ys[::-1], ys[1:]])
    X = np.empty(2 * len(xc) - 1); X[0::2] = xc; X[1::2] = 0.5 * (xc[:-1] + xc[1:])
    Y = np.empty(2 * len(yc) - 1); Y[0::2] = yc; Y[1::2] = 0.5 * (yc[:-1] + yc[1:])
    nx, ny = len(X), len(Y)
    j0 = 2 * (len(ys) - 1)
    iL = int(np.argmin(np.abs(X + 2 * a))); i0 = int(np.argmin(np.abs(X)))
    assert abs(X[iL] + 2 * a) < 1e-12 and abs(X[i0]) < 1e-12
    nid = np.arange(nx * ny).reshape(ny, nx)
    dup, nn = {}, nx * ny
    for i in range(iL + 1, i0):
        dup[i] = nn; nn += 1
    coords = np.zeros((nn, 2))
    coords[:nx * ny, 0] = np.tile(X, ny); coords[:nx * ny, 1] = np.repeat(Y, nx)
    for i, k in dup.items():
        coords[k] = (X[i], 0.0)
    ne_x, ne_y = len(xc) - 1, len(yc) - 1
    conn = np.zeros((ne_x * ne_y, 9), dtype=np.int64); mat = np.zeros(ne_x * ne_y, dtype=np.int64)
    e = 0
    for ey in range(ne_y):
        for ex in range(ne_x):
            nodes = []
            for b in range(3):
                for aa in range(3):
                    j, i = 2 * ey + b, 2 * ex + aa
                    k = nid[j, i]
                    if j == j0 and i in dup and ey < len(ys) - 1:
                        k = dup[i]
                    nodes.append(k)
            conn[e] = nodes; mat[e] = 1 if ey >= len(ys) - 1 else 2; e += 1
    hx = (xc[1:] - xc[:-1])[np.tile(np.arange(ne_x), ne_y)]
    hy = (yc[1:] - yc[:-1])[np.repeat(np.arange(ne_y), ne_x)]
    D = {1: plane_strain_D(E1, nu1), 2: plane_strain_D(E2, nu2)}
    Ke = np.zeros((len(conn), 18, 18))
    for tp, wp in zip(GP, GW):
        for tq, wq in zip(GP, GW):
            Na, Nb, dNa, dNb = n1d(tp), n1d(tq), dn1d(tp), dn1d(tq)
            dNdx = np.array([dNa[i] * Nb[j] for j in range(3) for i in range(3)])[None, :] * (2 / hx)[:, None]
            dNdy = np.array([Na[i] * dNb[j] for j in range(3) for i in range(3)])[None, :] * (2 / hy)[:, None]
            B = np.zeros((len(conn), 3, 18))
            B[:, 0, 0::2] = dNdx; B[:, 1, 1::2] = dNdy; B[:, 2, 0::2] = dNdy; B[:, 2, 1::2] = dNdx
            for m in (1, 2):
                s = mat == m
                Ke[s] += np.einsum('eki,kl,elj->eij', B[s], D[m], B[s]) * (wp * wq * hx[s] * hy[s] / 4)[:, None, None]
    dofs = np.zeros((len(conn), 18), dtype=np.int64); dofs[:, 0::2] = 2 * conn; dofs[:, 1::2] = 2 * conn + 1
    K = sp.csr_matrix((Ke.ravel(), (np.repeat(dofs, 18, axis=1).ravel(), np.tile(dofs, (1, 18)).ravel())), shape=(2 * nn, 2 * nn))
    # remote stress: sigma11 in material 2 zero; material 1 from continuity of eps11 (plane strain)
    s11_1 = (nu1 * (1 + nu1) / E1 - nu2 * (1 + nu2) / E2) * sig * E1 / (1 - nu1 ** 2)
    F = np.zeros(2 * nn)
    def edge(nodes3, p0, p1, t):          # consistent nodal forces of a uniform traction t on a 3-node edge
        L = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
        for tq, wq in zip(GP, GW):
            Nb = n1d(tq)
            for b in range(3):
                F[2 * nodes3[b]] += t[0] * Nb[b] * wq * L / 2
                F[2 * nodes3[b] + 1] += t[1] * Nb[b] * wq * L / 2
    for ex in range(ne_x):                 # top (n = e2) and bottom (n = -e2)
        edge([nid[ny - 1, 2 * ex + k] for k in range(3)], (xc[ex], W), (xc[ex + 1], W), (tau, sig))
        edge([nid[0, 2 * ex + k] for k in range(3)], (xc[ex], -W), (xc[ex + 1], -W), (-tau, -sig))
    for ey in range(ne_y):                 # right (n = e1) and left (n = -e1)
        s11 = s11_1 if ey >= len(ys) - 1 else 0.0
        edge([nid[2 * ey + k, nx - 1] for k in range(3)], (W, yc[ey]), (W, yc[ey + 1]), (s11, tau))
        edge([nid[2 * ey + k, 0] for k in range(3)], (-W, yc[ey]), (-W, yc[ey + 1]), (-s11, -tau))
    # rigid motion: fix the node at the bottom-left corner fully and the bottom-right corner vertically
    fixed = [2 * nid[0, 0], 2 * nid[0, 0] + 1, 2 * nid[0, nx - 1] + 1]
    free = np.setdiff1d(np.arange(2 * nn), fixed)
    u = np.zeros(2 * nn)
    u[free] = spla.spsolve(K[free][:, free].tocsc(), F[free], permc_spec='COLAMD')
    alpha, beta, eps, c1, c2 = dundurs(E1, nu1, E2, nu2)
    ii = np.arange(i0 - 1, iL, -1)
    r = -X[ii]
    up = np.array([nid[j0, i] for i in ii]); lo = np.array([dup[i] for i in ii])
    d1 = u[2 * up] - u[2 * lo]; d2 = u[2 * up + 1] - u[2 * lo + 1]
    fac = (c1 + c2) / (2 * math.sqrt(2 * math.pi) * (1 + 2j * eps) * math.cosh(math.pi * eps))
    Kapp = (d2 + 1j * d1) / (fac * np.sqrt(r) * r ** (1j * eps))
    sel = (r >= 4e-3 * a) & (r <= 4e-2 * a)
    A = np.vstack([np.ones(sel.sum()), r[sel], r[sel] ** 2]).T
    Kfe = np.linalg.lstsq(A, Kapp[sel], rcond=None)[0][0]
    Kex = (sig + 1j * tau) * (1 + 2j * eps) * math.sqrt(math.pi * a) * (2 * a) ** (-1j * eps)
    return alpha, beta, eps, Kfe, Kex, 2 * nn


if __name__ == '__main__':
    bad = 0
    cases = [(5 * 1.35 / 1.22 * 1.22 / 1.35, 0.22, 1.0, 0.35, 1.0, 0.0),   # E1/E2 = 5, nu1 0.22, nu2 0.35: tension
             (5.0, 0.22, 1.0, 0.35, 0.0, 1.0),                                # same materials, shear
             (1 / 6.0, 0.45, 1.0, 0.2, 1.0, 0.5)]                             # compliant on stiff, mixed
    for E1, nu1, E2, nu2, sg, ta in cases:
        al, be, ep, Kfe, Kex, nd = run(E1, nu1, E2, nu2, sg, ta)
        dpsi = math.degrees(np.angle(Kfe / Kex))
        dmag = abs(Kfe) / abs(Kex) - 1
        ok = abs(dpsi) < 0.1 and abs(dmag) < 0.01
        bad += not ok
        print('%s alpha %+.4f beta %+.4f eps %+.5f | K exact %.5f%+.5fi  FE %.5f%+.5fi | phase diff %+.4f deg, |K| diff %+.4f %% | ndof %d' %
              ('ok  ' if ok else 'FAIL', al, be, ep, Kex.real, Kex.imag, Kfe.real, Kfe.imag, dpsi, 100 * dmag, nd), flush=True)
    sys.exit(1 if bad else 0)
