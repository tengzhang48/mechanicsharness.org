"""
Independent FEM solver: semi-infinite interface crack in a bilayer (arbiter).

Plane strain, bilinear Q4 quads, structured mesh graded toward the crack tip
(x=0) in x and toward the interface (y=0) in y.  The interface nodes are
shared for x >= 0 (bonded) and duplicated for x < 0 (crack faces).

Loads: far-beam tractions sigma11(y) applied at x = -L as consistent nodal
loads (sigma12 = sigma22 = 0 there); all other external faces free.

K extraction: average stresses of the two elements adjacent to the interface
(evaluated at their centres), ahead of the tip; fit
   (sigma22 + i sigma12)(r) sqrt(2 pi r) e^{-i eps ln r} = K + c1 r + c2 r^2.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def make_mesh(h, H, L, dmin, growth, ny1, ny2, yfact):
    # --- x lines: graded from 0 outward
    right = [0.0]
    x, d = 0.0, dmin
    while x < L:
        x = min(x + d, L)
        right.append(x)
        d = min(d * growth, 0.45 * L)
    xs = np.array([-v for v in right[1:]][::-1] + right)

    # --- y lines above 0 (layer 1): graded from the interface
    up = [0.0]
    y, d = 0.0, dmin * yfact
    while y < h:
        y = min(y + d, h)
        up.append(y)
        d = min(d * growth, h / ny1)
    if len(up) - 1 < ny1:                      # ensure minimum layers
        up = list(np.linspace(0, h, ny1 + 1))
    dn = [0.0]
    y, d = 0.0, dmin * yfact
    while y < H:
        y = min(y + d, H)
        dn.append(y)
        d = min(d * growth, H / ny2)
    if len(dn) - 1 < ny2:
        dn = list(np.linspace(0, H, ny2 + 1))
    ys_up = np.array(up)
    ys_dn = -np.array(dn)[::-1]           # negative y, ascending: -H ... 0
    return xs, ys_dn, ys_up


def solve_fem(h, H, L, E1, E2, nu, P, M, dmin=2e-4, growth=1.28,
              ny1=16, ny2=24, yfact=0.5, fit_rmin=None, fit_rmax=None,
              nterms=3, return_fields=False):
    mu1, mu2 = E1 / (2 * (1 + nu)), E2 / (2 * (1 + nu))
    kap = 3 - 4 * nu
    Lam1, Lam2 = 2 * mu1 * nu / (1 - 2 * nu), 2 * mu2 * nu / (1 - 2 * nu)

    xs, ys_dn, ys_up = make_mesh(h, H, L, dmin, growth, ny1, ny2, yfact)
    nx = len(xs) - 1
    ndn = len(ys_dn) - 1        # element rows below interface
    nup = len(ys_up) - 1        # element rows above
    ix0 = int(np.searchsorted(xs, 0.0))          # first node line with x>=0

    # node ids: lower block rows iy=0..ndn (ndn+1 rows), upper block iy=0..nup
    ns = nx + 1                 # nodes per row
    def id_l(ix, iy):    # lower block (includes its interface row iy=ndn)
        return iy * ns + ix
    base_u = (ndn + 1) * ns
    def id_u(ix, iy):    # upper block (iy=0 is its interface row)
        return base_u + iy * ns + ix
    n_nodes_l = (ndn + 1) * ns
    n_nodes_u = (nup + 1) * ns
    # merge map: upper interface node ix>=ix0 -> lower interface node
    merge = {}
    for ix in range(ix0, nx):
        merge[id_u(ix, 0)] = id_l(ix, ndn)

    def dof(n):
        n = merge.get(n, n)
        return 2 * n

    # aliased upper-interface nodes (ix >= ix0) keep their ids but are never
    # referenced; size the system for all ids and pin the dead ones below.
    ndof = 2 * (n_nodes_l + n_nodes_u)

    D1 = np.array([[Lam1 + 2 * mu1, Lam1, 0], [Lam1, Lam1 + 2 * mu1, 0], [0, 0, mu1]])
    D2 = np.array([[Lam2 + 2 * mu2, Lam2, 0], [Lam2, Lam2 + 2 * mu2, 0], [0, 0, mu2]])

    gp = (-1 / np.sqrt(3), 1 / np.sqrt(3))
    rows, cols, vals = [], [], []

    def elem(ix, iy, ya, yb, upper):
        xa, xb = xs[ix], xs[ix + 1]
        D = D1 if upper else D2
        idf = id_u if upper else id_l
        iy_lo = iy if upper else iy          # upper: iy=0 interface row
        nodes = [idf(ix, iy), idf(ix + 1, iy), idf(ix + 1, iy + 1), idf(ix, iy + 1)]
        dofs = []
        for nd in nodes:
            ndm = merge.get(nd, nd)
            dofs += [2 * ndm, 2 * ndm + 1]
        ke = np.zeros((8, 8))
        for r in gp:
            for s in gp:
                dNr = np.array([-(1 - s), (1 - s), (1 + s), -(1 + s)]) / 4
                dNs = np.array([-(1 - r), -(1 + r), (1 + r), (1 - r)]) / 4
                detJ = (xb - xa) * (yb - ya) / 4
                dNx, dNy = dNr * 2 / (xb - xa), dNs * 2 / (yb - ya)
                B = np.zeros((3, 8))
                B[0, 0::2] = dNx
                B[1, 1::2] = dNy
                B[2, 0::2] = dNy
                B[2, 1::2] = dNx
                ke += B.T @ D @ B * detJ
        for a in range(8):
            for b in range(8):
                rows.append(dofs[a]); cols.append(dofs[b]); vals.append(ke[a, b])

    for iy in range(ndn):
        for ix in range(nx):
            elem(ix, iy, ys_dn[iy], ys_dn[iy + 1], upper=False)
    for iy in range(nup):
        for ix in range(nx):
            elem(ix, iy, ys_up[iy], ys_up[iy + 1], upper=True)

    Kmat = sp.csr_matrix((vals, (rows, cols)), shape=(ndof, ndof))

    # consistent loads at x = -L: traction t = sigma11 e_x on the end face
    # (outward normal -x, so applied nodal force = -int sigma11 N dA)
    f = np.zeros(ndof)

    def sig1(y, upper):
        if upper:
            return -P / h - 12 * M * (y - h / 2) / h ** 3
        return P / H + 12 * (M + P * (h + H) / 2) * (y + H / 2) / H ** 3

    def edge_load(ya, yb, n1, n2, upper):
        # 2-point Gauss along the edge for -int sigma11 N dy
        for g in gp:
            yy = ya + (g + 1) / 2 * (yb - ya)
            w = (yb - ya) / 2
            N1, N2 = (1 - g) / 2, (1 + g) / 2
            for (nd, N) in [(n1, N1), (n2, N2)]:
                ndm = merge.get(nd, nd)
                f[2 * ndm] += -sig1(yy, upper) * N * w

    for iy in range(ndn):
        edge_load(ys_dn[iy], ys_dn[iy + 1], id_l(0, iy), id_l(0, iy + 1), upper=False)
    for iy in range(nup):
        edge_load(ys_up[iy], ys_up[iy + 1], id_u(0, iy), id_u(0, iy + 1), upper=True)

    # pin rigid body and dead (aliased-away) dofs
    used = np.zeros(ndof, dtype=bool)
    used[np.array(rows, dtype=int)] = True
    dead = np.where(~used)[0]
    fix = list(dead) + [2 * id_l(0, 0), 2 * id_l(0, 0) + 1, 2 * id_l(nx - 1, 0) + 1]
    Kl = Kmat.tolil()
    for d in fix:
        Kl[d, :] = 0.0
        Kl[d, d] = 1.0
        f[d] = 0.0
    Kmat = Kl.tocsr()
    u = spla.spsolve(Kmat, f)

    # ---- K extraction ----------------------------------------------------
    eps = np.log((kap / mu1 + 1 / mu2) / (kap / mu2 + 1 / mu1)) / (2 * np.pi)
    rs, ts = [], []
    for ix in range(ix0, nx - 1):
        xa, xb = xs[ix], xs[ix + 1]
        xc = 0.5 * (xa + xb)
        ya_d = 0.5 * ys_dn[ndn - 1]           # centre of top lower element row? no:
        # centres of the elements adjacent to the interface:
        ycl = 0.5 * (ys_dn[ndn - 1] + ys_dn[ndn])      # top lower element centre
        ycu = 0.5 * (ys_up[0] + ys_up[1])              # bottom upper element centre
        # compute stresses at the two centres
        def stress_at(xc_, yc_, upper):
            iy = (ndn - 1) if not upper else 0
            ix_ = min(np.searchsorted(xs, xc_) - 1, nx - 2)
            ya, yb = (ys_dn[iy], ys_dn[iy + 1]) if not upper else (ys_up[iy], ys_up[iy + 1])
            xa, xb = xs[ix_], xs[ix_ + 1]
            idf = id_u if upper else id_l
            nodes = [idf(ix_, iy), idf(ix_ + 1, iy), idf(ix_ + 1, iy + 1), idf(ix_, iy + 1)]
            uv = []
            for nd in nodes:
                ndm = merge.get(nd, nd)
                uv += [u[2 * ndm], u[2 * ndm + 1]]
            uv = np.array(uv)
            r_, s_ = (2 * xc_ - xa - xb) / (xb - xa), (2 * yc_ - ya - yb) / (yb - ya)
            dNr = np.array([-(1 - s_), (1 - s_), (1 + s_), -(1 + s_)]) / 4
            dNs = np.array([-(1 - r_), -(1 + r_), (1 + r_), (1 - r_)]) / 4
            dNx, dNy = dNr * 2 / (xb - xa), dNs * 2 / (yb - ya)
            B = np.zeros((3, 8))
            B[0, 0::2] = dNx
            B[1, 1::2] = dNy
            B[2, 0::2] = dNy
            B[2, 1::2] = dNx
            strn = B @ uv
            D = D1 if upper else D2
            return D @ strn
        sl = stress_at(xc, ycl, upper=False)
        su = stress_at(xc, ycu, upper=True)
        sig22 = 0.5 * (sl[1] + su[1])
        sig12 = 0.5 * (sl[2] + su[2])
        rs.append(xc)
        ts.append(sig22 + 1j * sig12)
    rs = np.array(rs); ts = np.array(ts)
    if fit_rmin is None: fit_rmin = 4 * dmin
    if fit_rmax is None: fit_rmax = 0.25 * h
    m = (rs >= fit_rmin) & (rs <= fit_rmax)
    rfit, tfit = rs[m], ts[m]
    yv = tfit * np.sqrt(2 * np.pi * rfit) * np.exp(-1j * eps * np.log(rfit))
    Afit = np.vander(rfit, nterms, increasing=True)
    coef, *_ = np.linalg.lstsq(Afit, yv, rcond=None)
    K = coef[0]
    return dict(K=K, Kabs=abs(K), psi_deg=float(np.degrees(np.angle(K))),
                eps=float(eps), fit_r=(float(rfit[0]), float(rfit[-1])),
                nfit=len(rfit), ndof=ndof, nelem=2 * nx * (ndn + nup),
                resid=float(np.sqrt(np.mean(np.abs(yv - Afit @ coef) ** 2)) /
                           np.sqrt(np.mean(np.abs(yv) ** 2))),
                u=u, xs=xs, merge=merge, id_l=id_l, id_u=id_u,
                ys_dn=ys_dn, ys_up=ys_up)
