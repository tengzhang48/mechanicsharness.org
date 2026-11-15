"""Triangular mesh of a rectangle [X0,X1]x[Y0,Y1] with an interface on y = 0 and crack segments on y = 0.

Around each crack tip: a square hole of half-size R0 filled with a polar ("spider-web") mesh, rings
geometrically graded toward the tip (radial step ~ r * 2pi/Ntheta), Ntheta = 8*nh sectors, rings blended
from circles (inner) to the square (outer ring).  Outside the holes: a graded tensor grid of right
triangles (x and y node sets contain the uniform hole-side nodes, spacing d0 = R0/nh, and grow away from
the holes with ratio controlled by beta, capped by smax).  Crack faces carry duplicated nodes.
No 9-node quadrilaterals anywhere; P2 triangles are used by the solver on this P1 geometry.
"""
import numpy as np


def graded_nodes(A, B, size_fun, nsample=20001):
    """Nodes from A to B with local spacing ~ size_fun(x); exact end points."""
    if B <= A:
        return np.array([A])
    xs = np.linspace(A, B, nsample)
    inv = 1.0 / size_fun(xs)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(xs))])
    n = max(1, int(np.ceil(cum[-1] - 1e-9)))
    targets = np.linspace(0, cum[-1], n + 1)
    nodes = np.interp(targets, cum, xs)
    nodes[0], nodes[-1] = A, B
    return nodes


def build_mesh(X0, X1, Y0, Y1, tips, cracks, R0, nh, beta, smax_x, smax_y, rmin, rho_blend=0.5):
    """tips: list of (xt, side) with side=-1 if the crack lies at x<xt, +1 if at x>xt (tips on y=0).
    cracks: list of (xa, xb) crack intervals on y=0 (xa<xb), tips at their interior ends.
    Returns (points[N,2], cells[M,3], cell_mat[M] (1: y>0, 2: y<0))."""
    d0 = R0 / nh
    holes = sorted([(xt - R0, xt + R0) for xt, _ in tips])

    def dist_hole_x(x):
        d = np.full_like(x, np.inf)
        for a, b in holes:
            d = np.minimum(d, np.where(x < a, a - x, np.where(x > b, x - b, 0.0)))
        return d

    sx = lambda x: np.minimum(smax_x, d0 + beta * dist_hole_x(x))
    # x nodes
    fixed = [X0]
    for a, b in holes:
        fixed += [a, b]
    fixed += [X1]
    xn = [np.array([X0])]
    for k in range(len(fixed) - 1):
        A, B = fixed[k], fixed[k + 1]
        if k % 2 == 1:  # inside a hole: uniform
            seg = np.linspace(A, B, 2 * nh + 1)
        else:
            seg = graded_nodes(A, B, sx)
        xn.append(seg[1:])
    xn = np.concatenate(xn)
    # y nodes
    sy = lambda y: np.minimum(smax_y, d0 + beta * np.maximum(np.abs(y) - R0, 0.0))
    ylo = graded_nodes(Y0, -R0, sy)
    ymid = np.linspace(-R0, R0, 2 * nh + 1)
    yhi = graded_nodes(R0, Y1, sy)
    yn = np.concatenate([ylo, ymid[1:], yhi[1:]])
    yn[np.argmin(np.abs(yn))] = 0.0
    j0 = int(np.argmin(np.abs(yn)))

    def in_crack(x):
        return any(a + 1e-12 < x < b - 1e-12 for a, b in cracks)

    def in_hole(xc, yc):
        return any(a < xc < b for a, b in holes) and (-R0 < yc < R0)

    nodes = {}
    pts = []

    def nid(key, xy):
        if key not in nodes:
            nodes[key] = len(pts)
            pts.append(xy)
        return nodes[key]

    def tensor_node(i, j, upper):
        flag = 0
        if j == j0 and in_crack(xn[i]):
            flag = 1 if upper else -1
        return nid(("t", i, j, flag), (xn[i], yn[j]))

    cells = []
    for i in range(len(xn) - 1):
        for j in range(len(yn) - 1):
            xc, yc = 0.5 * (xn[i] + xn[i + 1]), 0.5 * (yn[j] + yn[j + 1])
            if in_hole(xc, yc):
                continue
            up = yc > 0
            n00, n10 = tensor_node(i, j, up), tensor_node(i + 1, j, up)
            n11, n01 = tensor_node(i + 1, j + 1, up), tensor_node(i, j + 1, up)
            cells.append((n00, n10, n11))
            cells.append((n00, n11, n01))

    # polar zones
    Nt = 8 * nh
    q = np.exp(2 * np.pi / Nt)
    Nr = int(np.ceil(np.log(R0 / rmin) / np.log(q)))
    rho = q ** (np.arange(1, Nr + 1) - Nr)  # rho[-1] = 1
    # square boundary in local frame (crack along local -x), counterclockwise from (-R0, 0-)
    S = []
    for k in range(nh):
        S.append((-R0, -k * d0))
    for k in range(2 * nh):
        S.append((-R0 + k * d0, -R0))
    for k in range(2 * nh):
        S.append((R0, -R0 + k * d0))
    for k in range(2 * nh):
        S.append((R0 - k * d0, R0))
    for k in range(nh + 1):
        S.append((-R0, R0 - k * d0))
    S = np.array(S)
    assert len(S) == Nt + 1
    thu = -np.pi + 2 * np.pi * np.arange(Nt + 1) / Nt
    ix_of = lambda x: int(np.argmin(np.abs(xn - x)))
    iy_of = lambda y: int(np.argmin(np.abs(yn - y)))
    for (xt, side) in tips:
        sgn = 1.0 if side == -1 else -1.0  # rotation by pi when the crack lies at x > xt
        tip = nid(("tip", xt), (xt, 0.0))
        P = np.zeros((Nr + 1, Nt + 1), dtype=int)
        for ir in range(1, Nr + 1):
            rr = rho[ir - 1]
            w = 0.0 if rr <= rho_blend else ((rr - rho_blend) / (1 - rho_blend)) ** 2
            for jt in range(Nt + 1):
                loc = rr * ((1 - w) * R0 * np.array([np.cos(thu[jt]), np.sin(thu[jt])]) + w * S[jt])
                if jt in (0, Nt, Nt // 2):
                    loc[1] = 0.0
                g = np.array([xt + sgn * loc[0], sgn * loc[1]])
                # local lower face (jt=0) -> global flag: -1 if side=-1, +1 if rotated
                if jt == 0:
                    flag = -1 if side == -1 else 1
                elif jt == Nt:
                    flag = 1 if side == -1 else -1
                else:
                    flag = 0
                if ir == Nr:
                    gx, gy = xt + sgn * S[jt][0], sgn * S[jt][1]
                    i, j = ix_of(gx), iy_of(gy)
                    assert abs(xn[i] - gx) < 1e-9 * max(1, abs(gx)) and abs(yn[j] - gy) < 1e-9 * max(1, abs(gy))
                    upper = (flag == 1)
                    if flag == 0:
                        P[ir, jt] = nid(("t", i, j, 0), (xn[i], yn[j]))
                    else:
                        P[ir, jt] = tensor_node(i, j, upper)
                else:
                    if flag == 0:
                        P[ir, jt] = nid(("p", xt, ir, jt), tuple(g))
                    else:
                        P[ir, jt] = nid(("p", xt, ir, jt, flag), tuple(g))
        for jt in range(Nt):
            cells.append((tip, P[1, jt], P[1, jt + 1]))
        for ir in range(1, Nr):
            for jt in range(Nt):
                a, b, c, d = P[ir, jt], P[ir + 1, jt], P[ir + 1, jt + 1], P[ir, jt + 1]
                cells.append((a, b, c))
                cells.append((a, c, d))
    pts = np.array(pts, float)
    cells = np.array(cells, dtype=np.int64)
    cen = pts[cells].mean(axis=1)
    mat = np.where(cen[:, 1] > 0, 1, 2)
    info = dict(xn=xn, yn=yn, Nr=Nr, Nt=Nt, q=q, d0=d0)
    return pts, cells, mat, info


def check_mesh(pts, cells, X0, X1, Y0, Y1, crack_len):
    p = pts[cells]
    area = 0.5 * ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                  - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    edges = np.sort(np.concatenate([cells[:, [0, 1]], cells[:, [1, 2]], cells[:, [2, 0]]]), axis=1)
    e, cnt = np.unique(edges, axis=0, return_counts=True)
    be = e[cnt == 1]
    blen = np.sum(np.linalg.norm(pts[be[:, 0]] - pts[be[:, 1]], axis=1))
    expect_area = (X1 - X0) * (Y1 - Y0)
    expect_blen = 2 * (X1 - X0) + 2 * (Y1 - Y0) + 2 * crack_len
    used = np.zeros(len(pts), bool)
    used[cells.ravel()] = True
    # aspect ratio (longest edge^2 / (2*area))
    L = np.stack([np.linalg.norm(p[:, 1] - p[:, 0], axis=1), np.linalg.norm(p[:, 2] - p[:, 1], axis=1),
                  np.linalg.norm(p[:, 0] - p[:, 2], axis=1)], axis=1)
    asp = L.max(axis=1) ** 2 / (2 * area)
    return dict(min_area=area.min(), area_err=abs(area.sum() - expect_area) / expect_area,
                blen_err=abs(blen - expect_blen) / expect_blen, max_edge_count=cnt.max(),
                unused_nodes=int((~used).sum()), ncells=len(cells), npts=len(pts), max_aspect=asp.max())
