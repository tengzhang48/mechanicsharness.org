"""Finite-element calculation for the semi-infinite bilayer interface crack.

The only external unit choices are E2=h=1 and H=2.  The crack tip is (0,0).
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

NU = 1.0 / 3.0
KAPPA = 3.0 - 4.0 * NU
H = 2.0
E1, E2 = 9.0, 1.0
MU1 = E1 / (2.0 * (1.0 + NU))
MU2 = E2 / (2.0 * (1.0 + NU))
EPS = math.log((KAPPA / MU1 + 1.0 / MU2) /
               (KAPPA / MU2 + 1.0 / MU1)) / (2.0 * math.pi)
LAMBDA = 0.5 + 1j * EPS
S = (KAPPA + 1.0) / MU1 + (KAPPA + 1.0) / MU2
G_FACTOR = S / (16.0 * math.cosh(math.pi * EPS) ** 2)  # G = G_FACTOR |K|^2
OPEN_FACTOR = S / (2.0 * math.sqrt(2.0 * math.pi) *
                   math.cosh(math.pi * EPS) * (1.0 + 2j * EPS))


def graded_axis(length, first, growth, largest):
    """Nodes 0..length, geometrically growing away from the first node."""
    steps = [first]
    while True:
        end = sum(steps)
        if end >= length * (1.0 - 1e-13):
            break
        nxt = min(steps[-1] * growth, largest)
        if end + nxt > length:
            nxt = length - end
        steps.append(nxt)
        if len(steps) > 10000:
            raise RuntimeError('axis generation failed')
    y = np.r_[0.0, np.cumsum(steps)]
    y[-1] = length
    return y


def make_grid(length, first, growth, largest_x, largest_y):
    xp = graded_axis(length, first, growth, largest_x)
    x = np.r_[-xp[:0:-1], xp]
    yt_int = graded_axis(1.0, first, growth, largest_y)
    yb_int = graded_axis(H, first, growth, largest_y)
    yt = yt_int
    yb = -yb_int[::-1]
    return x, yb, yt


def assemble_material(x, y, bottom, E):
    """Assemble one rectangular layer.

    bottom=True means interface is the upper grid line, otherwise lower.
    Returns (element connectivity, stiffness matrices).  Node numbers are
    local row-major indices with x varying fastest.
    """
    nx, ny = len(x), len(y)
    ix, iy = np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), indexing='xy')
    ix, iy = ix.ravel(), iy.ravel()
    n0 = iy * nx + ix
    conn = np.column_stack((n0, n0 + 1, n0 + 1 + nx, n0 + nx))
    hx = x[ix + 1] - x[ix]
    hy = y[iy + 1] - y[iy]
    D = E / ((1 + NU) * (1 - 2 * NU)) * np.array(
        [[1 - NU, NU, 0], [NU, 1 - NU, 0], [0, 0, (1 - 2 * NU) / 2]])
    ke = np.zeros((len(conn), 8, 8))
    gp = 1.0 / math.sqrt(3.0)
    for xi in (-gp, gp):
        for eta in (-gp, gp):
            dN_dxi = 0.25 * np.array(
                [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
            dN_deta = 0.25 * np.array(
                [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
            dNdx = dN_dxi[None, :] * (2.0 / hx[:, None])
            dNdy = dN_deta[None, :] * (2.0 / hy[:, None])
            B = np.zeros((len(conn), 3, 8))
            B[:, 0, 0::2] = dNdx
            B[:, 1, 1::2] = dNdy
            B[:, 2, 0::2] = dNdy
            B[:, 2, 1::2] = dNdx
            detJ = hx * hy / 4.0
            ke += np.einsum('eja,jk,ekb,e->eab', B, D, B, detJ)
    return conn, ke


def solve_model(length=12.0, first=1e-4, growth=1.16,
                largest_x=0.20, largest_y=0.15, verbose=True):
    t0 = time.time()
    x, yb, yt = make_grid(length, first, growth, largest_x, largest_y)
    nx, nb, nt = len(x), len(yb), len(yt)
    if verbose:
        print(f'grid: nx={nx}, bottom={nb}, top={nt}, '
              f'min dx={np.min(np.diff(x)):.3g}', flush=True)

    cb, kb = assemble_material(x, yb, True, E2)
    ct, kt = assemble_material(x, yt, False, E1)
    # Bottom nodes are kept.  Top nodes are appended, except bonded interface
    # nodes at x>=0, which map to the existing upper line of bottom nodes.
    top_map = np.arange(nx * nb, nx * nb + nx * nt, dtype=np.int64)
    tip = nx // 2
    top_map[:nx] = np.arange(nx * nb, nx * nb + nx)
    top_map[tip:nx] = (nb - 1) * nx + np.arange(tip, nx)
    nnode = nx * nb + nx * nt
    ndof = 2 * nnode

    # Account for unused duplicate top-interface nodes on the bonded side by
    # simply retaining them with zero stiffness, then eliminating them below.
    conn_all = np.vstack((cb, top_map[ct]))
    ke_all = np.vstack((kb, kt.reshape(-1, 8, 8))) if kt.ndim != 3 else np.vstack((kb, kt))
    edofs = np.empty((len(conn_all), 8), dtype=np.int64)
    edofs[:, 0::2] = 2 * conn_all
    edofs[:, 1::2] = 2 * conn_all + 1
    rows = np.repeat(edofs, 8, axis=1).ravel()
    cols = np.tile(edofs, (1, 8)).ravel()
    K = sp.coo_matrix((ke_all.ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()
    K.sum_duplicates()

    # Consistent tractions on the far left ends of the detached layers.
    # Two load vectors: case A (P=1,M=0), case B (P=0,M=1).
    F = np.zeros((ndof, 2))
    gp = 1 / math.sqrt(3)
    for y, E in ((yb, E2), (yt, E1)):
        for j in range(len(y) - 1):
            y0, y1 = y[j:j + 2]
            for eta in (-gp, gp):
                ym = 0.5 * ((1 - eta) * y0 + (1 + eta) * y1)
                jac = 0.5 * (y1 - y0)
                N = np.array([(1 - eta) / 2, (1 + eta) / 2])
                if y is yt:
                    sA = -1.0
                    sB = -12.0 * (ym - 0.5)
                else:
                    sA = 1.0 / H + 12.0 * ((H + 1) / 2) * (ym + H / 2) / H**3
                    sB = 12.0 * (ym + H / 2) / H**3
                nodes = np.array([j * nx, (j + 1) * nx])
                if y is yt:
                    nodes = top_map[np.array([j * nx, (j + 1) * nx])]
                F[2 * nodes, 0] -= N * sA * jac
                F[2 * nodes, 1] -= N * sB * jac

    # Isolated duplicate bonded top-interface nodes were never used.
    used = np.zeros(nnode, dtype=bool)
    used[conn_all.ravel()] = True
    # Minimal rigid-body restraints at the far bonded end.
    interface_right = (nb - 1) * nx + nx - 1
    top_right = top_map[(nt - 1) * nx + nx - 1]
    used[interface_right] = True
    used[top_right] = True
    restrained = {2 * interface_right, 2 * interface_right + 1,
                  2 * top_right + 1}
    active_nodes = np.flatnonzero(used)
    active_dofs = np.r_[2 * active_nodes, 2 * active_nodes + 1]
    active_dofs.sort()
    free = np.array([d for d in active_dofs if d not in restrained], dtype=np.int64)
    Kr = K[free][:, free].tocsc()
    Fr = F[free]
    if verbose:
        print(f'unknowns={len(free)}, nnz={Kr.nnz}, assembly={time.time()-t0:.1f}s', flush=True)
        print('factorizing...', flush=True)
    lu = spla.splu(Kr, permc_spec='COLAMD')
    Ur = lu.solve(Fr)
    U = np.zeros((ndof, 2))
    U[free] = Ur
    if verbose:
        print(f'solved in {time.time()-t0:.1f}s', flush=True)

    # Crack-opening estimates at every negative interface grid node.
    bot_int = (nb - 1) * nx + np.arange(tip)
    top_int = top_map[np.arange(tip)]
    ub = U[2 * bot_int]
    vb = U[2 * bot_int + 1]
    ut = U[2 * top_int]
    vt = U[2 * top_int + 1]
    jump = (vt - vb) + 1j * (ut - ub)
    r = -x[:tip]
    Kline = jump / (OPEN_FACTOR * r[:, None] ** LAMBDA)
    return dict(x=x, yb=yb, yt=yt, r=r, Kline=Kline, U=U, F=F,
                top_map=top_map, nx=nx, nb=nb, nt=nt, free=free,
                eps=EPS, G_FACTOR=G_FACTOR)


def extrapolate_openings(result, rmax=0.08, rmin=0.0, degree=3):
    r = result['r']
    kl = result['Kline']
    mask = (r <= rmax) & (r > max(rmin, 0))
    # Exclude the singular-tip-adjacent point only if requested.
    rr = r[mask]
    z = np.sqrt(rr / rmax)
    out = []
    for j in range(kl.shape[1]):
        A = np.vander(z, degree + 1, increasing=True)
        c = np.linalg.lstsq(A, kl[mask, j], rcond=None)[0]
        out.append(c[0])
    return np.array(out)


def analytic_energy():
    # Plane-strain beam modulus E/(1-nu^2).
    eb1, eb2 = E1 / (1 - NU**2), E2 / (1 - NU**2)
    GA = 1 / (2 * eb1) + (1 / H + 12 * ((H + 1) / 2)**2 / H**3) / (2 * eb2)
    GB = 12 / (2 * eb1) + 12 / H**3 / (2 * eb2)
    return GA, GB


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--length', type=float, default=12.0)
    ap.add_argument('--first', type=float, default=1e-4)
    ap.add_argument('--growth', type=float, default=1.16)
    ap.add_argument('--largest-x', type=float, default=0.20)
    ap.add_argument('--largest-y', type=float, default=0.15)
    ap.add_argument('--out', type=str, default='')
    args = ap.parse_args()
    res = solve_model(args.length, args.first, args.growth,
                      args.largest_x, args.largest_y)
    GA, GB = analytic_energy()
    print(f'eps={EPS:.12g}, G factor={G_FACTOR:.12g}, open factor={OPEN_FACTOR}')
    print('analytic G, |K|, |K|sqrt terms:')
    for G in (GA, GB):
        Kabs = math.sqrt(G / G_FACTOR)
        print(f'G={G:.12g}, |K|={Kabs:.12g}')
    for rmax in (0.01, 0.02, 0.04, 0.08, 0.16, 0.32):
        for degree in (1, 2, 3):
            k = extrapolate_openings(res, rmax=rmax, degree=degree)
            print(f'rmax={rmax:g}, degree={degree}:',
                  ['|K|=%.9g psi=%.9g K=(%.9g,%.9g)' %
                   (abs(z), np.degrees(np.angle(z)), z.real, z.imag) for z in k])
    if args.out:
        # Save only small extraction data, not the large displacement vector.
        np.savez(args.out, r=res['r'], Kline=res['Kline'], x=res['x'],
                 yb=res['yb'], yt=res['yt'])


if __name__ == '__main__':
    main()
