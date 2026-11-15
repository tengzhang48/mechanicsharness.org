#!/usr/bin/env python3
"""Plane-strain FE calculation for the cracked finite-thickness bilayer.

The outer radius of every polar ray is the intersection with the rectangular
computational boundary.  This permits geometric refinement towards the crack
tip without a separate far-field mesh.  The two crack faces are independent;
interface nodes ahead of the tip are shared.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class Mesh:
    xy: np.ndarray
    top_elements: np.ndarray
    bottom_elements: np.ndarray
    top_triangles: np.ndarray
    bottom_triangles: np.ndarray
    top_outer: np.ndarray
    bottom_outer: np.ndarray
    crack_top: np.ndarray
    crack_bottom: np.ndarray
    crack_r: np.ndarray
    right_tip: int
    right_top_corner: int


class PolarLayer:
    def __init__(self, nrad: int, nang: int, rmin: float, length: float,
                 thickness: float, upper: bool):
        self.nrad = nrad
        self.upper = upper
        if upper:
            angles = np.linspace(0.0, np.pi, nang + 1)
            corners = [math.atan2(thickness, length),
                       np.pi - math.atan2(thickness, length)]
        else:
            angles = np.linspace(np.pi, 2.0 * np.pi, nang + 1)
            corners = [np.pi + math.atan2(thickness, length),
                       2.0 * np.pi - math.atan2(thickness, length)]
        self.theta = np.unique(np.r_[angles, corners])
        self.nray = len(self.theta)
        c, s = np.cos(self.theta), np.sin(self.theta)
        outer = np.full(self.nray, np.inf)
        def bounded_distance(numerator, denominator, mask):
            distance = np.full_like(denominator, np.inf, dtype=float)
            np.divide(numerator, denominator, out=distance, where=mask)
            return distance
        if upper:
            outer = np.minimum(outer, bounded_distance(
                thickness, s, s > 1e-14))
        else:
            outer = np.minimum(outer, bounded_distance(
                -thickness, s, s < -1e-14))
        outer = np.minimum(outer, bounded_distance(
            length, c, c > 1e-14))
        outer = np.minimum(outer, bounded_distance(
            -length, c, c < -1e-14))
        self.outer_radius = outer
        # Geometric radial grading, separately normalized on each ray.
        a = np.linspace(0.0, 1.0, nrad)
        self.radius = rmin * (outer[:, None] / rmin) ** a[None, :]
        local_xy = np.empty((self.nray, nrad, 2))
        local_xy[:, :, 0] = self.radius * c[:, None]
        local_xy[:, :, 1] = self.radius * s[:, None]
        self.local_xy = local_xy
        # Local node zero is the tip. Ring nodes then follow.
        self.nodes = np.arange(self.nray * nrad).reshape(
            self.nray, nrad) + 1
        self.tip = 0
        self.outer = self.nodes[:, -1]
        quads = []
        for k in range(self.nray - 1):
            for j in range(nrad - 1):
                quads.append([self.nodes[k, j], self.nodes[k, j + 1],
                              self.nodes[k + 1, j + 1],
                              self.nodes[k + 1, j]])
        self.quads = np.asarray(quads, dtype=np.int64)
        tris = []
        for k in range(self.nray - 1):
            tris.append([self.tip, self.nodes[k, 0],
                         self.nodes[k + 1, 0]])
        self.triangles = np.asarray(tris, dtype=np.int64)
        self.local_xy_all = np.vstack((np.zeros((1, 2)),
                                       local_xy.reshape(-1, 2)))


def build_mesh(nrad=180, nang=48, rmin=1e-4, length=20.0,
               h=1.0, H=2.0) -> Mesh:
    top = PolarLayer(nrad, nang, rmin, length, h, True)
    bottom = PolarLayer(nrad, nang, rmin, length, H, False)
    nt = len(top.local_xy_all)
    nb = len(bottom.local_xy_all)
    xy = np.vstack((top.local_xy_all, bottom.local_xy_all))
    top_global = np.arange(nt)
    bottom_global = np.arange(nt, nt + nb)

    # Merge the tip and every node on the bonded ray theta=0/2pi.
    # Top ray zero and bottom ray nray-1 have identical radii and coordinates.
    bottom_global[bottom.tip] = top_global[top.tip]
    for j in range(nrad):
        bt = bottom.nodes[-1, j]
        tt = top.nodes[0, j]
        assert np.allclose(bottom.local_xy_all[bt], top.local_xy_all[tt],
                           rtol=0.0, atol=1e-12)
        bottom_global[bt] = top_global[tt]

    # Compact node numbers after merging.
    all_old = np.r_[top_global, bottom_global]
    unique, inverse = np.unique(all_old, return_inverse=True)
    top_map = inverse[:nt]
    bottom_map = inverse[nt:]
    xy = xy[unique]
    tq = top_map[top.quads]
    bq = bottom_map[bottom.quads]
    tt = top_map[top.triangles]
    bt = bottom_map[bottom.triangles]
    top_outer = top_map[top.outer]
    bottom_outer = bottom_map[bottom.outer]
    crack_top = top_map[top.nodes[-1, :]]
    crack_bottom = bottom_map[bottom.nodes[0, :]]
    crack_r = top.radius[-1, :].copy()
    assert np.allclose(crack_r, bottom.radius[0, :])
    right_tip = top_map[top.nodes[0, -1]]
    # The right top corner is the outer node closest to (length,h).
    d = np.sum((xy[top_outer] - np.array([length, h])) ** 2, axis=1)
    right_top_corner = top_outer[np.argmin(d)]
    assert np.allclose(xy[right_top_corner], [length, h], atol=1e-12)
    return Mesh(xy, tq, bq, tt, bt, top_outer, bottom_outer,
                crack_top, crack_bottom, crack_r,
                int(right_tip), int(right_top_corner))


def plane_strain_D(E: float, nu: float) -> np.ndarray:
    return E / ((1.0 + nu) * (1.0 - 2.0 * nu)) * np.array(
        [[1.0 - nu, nu, 0.0], [nu, 1.0 - nu, 0.0],
         [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0]])


def quad_stiffness(xy: np.ndarray, elements: np.ndarray,
                   D: np.ndarray) -> np.ndarray:
    """Element matrices for bilinear quads with 2x2 Gauss integration."""
    gp = 1.0 / math.sqrt(3.0)
    gauss = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)]
    pxy = xy[elements]
    ke = np.zeros((len(elements), 8, 8))
    for xi, eta in gauss:
        dN = 0.25 * np.array([
            [-(1.0 - eta), -(1.0 - xi)],
            [ (1.0 - eta), -(1.0 + xi)],
            [ (1.0 + eta),  (1.0 + xi)],
            [-(1.0 + eta),  (1.0 - xi)]])
        # jac[natural, physical] = d(x_physical)/d(xi_natural)
        jac = np.einsum('gi,egj->eij', dN, pxy)
        det = np.linalg.det(jac)
        if np.any(det <= 0.0):
            raise ValueError(f'Nonpositive quad Jacobian: {det.min()}')
        dxy = np.einsum('eij,gj->egi', np.linalg.inv(jac), dN)
        B = np.zeros((len(elements), 3, 8))
        B[:, 0, 0::2] = dxy[:, :, 0]
        B[:, 1, 1::2] = dxy[:, :, 1]
        B[:, 2, 0::2] = dxy[:, :, 1]
        B[:, 2, 1::2] = dxy[:, :, 0]
        ke += np.einsum('eki,kl,elj,e->eij', B, D, B, det)
    return ke


def triangle_stiffness(xy: np.ndarray, elements: np.ndarray,
                       D: np.ndarray) -> np.ndarray:
    p = xy[elements]
    x1, y1 = p[:, 0, 0], p[:, 0, 1]
    x2, y2 = p[:, 1, 0], p[:, 1, 1]
    x3, y3 = p[:, 2, 0], p[:, 2, 1]
    twoarea = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if np.any(twoarea <= 0):
        raise ValueError('Nonpositive triangle area')
    area = twoarea / 2.0
    B = np.zeros((len(elements), 3, 8))[:, :, :6]
    B[:, 0, 0::2] = np.stack((y2 - y3, y3 - y1, y1 - y2), axis=1)
    B[:, 1, 1::2] = np.stack((x3 - x2, x1 - x3, x2 - x1), axis=1)
    B[:, 2, 0::2] = B[:, 1, 1::2]
    B[:, 2, 1::2] = B[:, 0, 0::2]
    B /= twoarea[:, None, None]
    k6 = np.einsum('eki,kl,elj,e->eij', B, D, B, area)
    ke = np.zeros((len(elements), 8, 8))
    inds = np.arange(6)
    ke[:, inds[:, None], inds[None, :]] = k6
    return ke


def assemble_stiffness(mesh: Mesh, h=1.0, H=2.0, E1=9.0, E2=1.0,
                       nu=1.0 / 3.0):
    n = len(mesh.xy)
    blocks = []
    for els, E in [(mesh.top_elements, E1),
                   (mesh.bottom_elements, E2)]:
        blocks.append((els, quad_stiffness(mesh.xy, els,
                                           plane_strain_D(E, nu))))
    for els, E in [(mesh.top_triangles, E1),
                   (mesh.bottom_triangles, E2)]:
        ke = triangle_stiffness(mesh.xy, els, plane_strain_D(E, nu))
        blocks.append((els[:, :3], ke[:, :6, :6]))
    rows, cols, vals = [], [], []
    for els, ke in blocks:
        edof = np.empty((len(els), els.shape[1] * 2), dtype=np.int64)
        edof[:, 0::2] = 2 * els
        edof[:, 1::2] = 2 * els + 1
        rows.append(np.repeat(edof, edof.shape[1], axis=1).ravel())
        cols.append(np.tile(edof, (1, edof.shape[1])).ravel())
        vals.append(ke.ravel())
    K = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(2 * n, 2 * n)).tocsr()
    return K


def add_left_load(f: np.ndarray, mesh: Mesh, P: float, M: float,
                  h=1.0, H=2.0):
    """Consistent nodal loads for the specified far-behind beam stresses."""
    gp = math.sqrt(3.0 / 5.0)
    gauss = [(-gp, 5.0 / 9.0), (0.0, 8.0 / 9.0), (gp, 5.0 / 9.0)]
    for outer, upper in [(mesh.top_outer, True),
                         (mesh.bottom_outer, False)]:
        for a, b in zip(outer[:-1], outer[1:]):
            xa, ya = mesh.xy[a]
            xb, yb = mesh.xy[b]
            # Left boundary is vertical; ignore top, bottom and right edges.
            if abs(xa + length_global) > 1e-9 or abs(xb + length_global) > 1e-9:
                continue
            if abs(xa - xb) > 1e-9:
                continue
            assert (ya >= -1e-10 and yb <= h + 1e-10) if upper else \
                   (ya <= 1e-10 and yb >= -H - 1e-10)
            for xi, w in gauss:
                s = (xi + 1.0) / 2.0
                y = ya * (1.0 - s) + yb * s
                if upper:
                    stress = -P / h - 12.0 * M * (y - h / 2.0) / h ** 3
                else:
                    M2 = M + P * (h + H) / 2.0
                    stress = P / H + 12.0 * M2 * (y + H / 2.0) / H ** 3
                # Outward normal is (-1,0), so tx=-sigma_11.
                traction_x = -stress
                jac = abs(yb - ya) / 2.0
                for node, Ni in [(a, 1.0 - s), (b, s)]:
                    f[2 * node] += Ni * traction_x * w * jac


def solve_load(K, mesh: Mesh, P: float, M: float, h=1.0, H=2.0):
    n = len(mesh.xy)
    f = np.zeros(2 * n)
    add_left_load(f, mesh, P, M, h, H)
    # Three zero-resultant rigid-body constraints. In exact arithmetic their
    # reactions vanish because the applied load has zero force and moment.
    fixed = np.array([2 * mesh.right_tip, 2 * mesh.right_tip + 1,
                      2 * mesh.right_top_corner])
    free = np.setdiff1d(np.arange(2 * n), fixed, assume_unique=False)
    u = np.zeros(2 * n)
    u[free] = spla.spsolve(K[free][:, free], f[free])
    return u, f


def hoop_opening(mesh: Mesh, u: np.ndarray):
    ut = u.reshape(-1, 2)[mesh.crack_top]
    ub = u.reshape(-1, 2)[mesh.crack_bottom]
    delta = ut - ub
    return delta[:, 1] + 1j * delta[:, 0]


def run(nrad=180, nang=48, rmin=1e-4, length=20.0,
        out='fem_raw.json'):
    global length_global
    length_global = length
    start = time.time()
    mesh = build_mesh(nrad, nang, rmin, length)
    print(f'mesh: {len(mesh.xy)} nodes, {len(mesh.top_elements)+len(mesh.bottom_elements)} quads')
    Kmat = assemble_stiffness(mesh)
    print(f'stiffness assembled at {time.time()-start:.1f}s, nnz={Kmat.nnz}')
    nu = 1.0 / 3.0
    kappa = 3.0 - 4.0 * nu
    mu1 = 9.0 / (2.0 * (1.0 + nu))
    mu2 = 1.0 / (2.0 * (1.0 + nu))
    beta = (mu1 * (kappa - 1.0) - mu2 * (kappa - 1.0)) / \
           (mu1 * (kappa + 1.0) + mu2 * (kappa + 1.0))
    eps = math.log((1.0 - beta) / (1.0 + beta)) / (2.0 * math.pi)
    D = (kappa + 1.0) / mu1 + (kappa + 1.0) / mu2
    # Solve P=1 and M=1 separately; linearity then gives every combination.
    uA, fA = solve_load(Kmat, mesh, 1.0, 0.0)
    print(f'case A solved at {time.time()-start:.1f}s')
    uB, fB = solve_load(Kmat, mesh, 0.0, 1.0)
    print(f'case B solved at {time.time()-start:.1f}s')
    r = mesh.crack_r
    deltaA = hoop_opening(mesh, uA)
    deltaB = hoop_opening(mesh, uB)
    # Both exponent signs are saved so that extraction can independently check
    # the sign convention. Standard convention: T ~ K r^(-1/2+i epsilon).
    result = {'parameters': {'nrad': nrad, 'nang': nang, 'rmin': rmin,
                             'length': length, 'nodes': len(mesh.xy),
                             'beta': beta, 'epsilon': eps, 'D': D},
              'r': r.tolist(), 'deltaA_real': deltaA.real.tolist(),
              'deltaA_imag': deltaA.imag.tolist(),
              'deltaB_real': deltaB.real.tolist(),
              'deltaB_imag': deltaB.imag.tolist()}
    with open(out, 'w') as f:
        json.dump(result, f)
    print(f'wrote {out}; elapsed {time.time()-start:.1f}s')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nrad', type=int, default=180)
    parser.add_argument('--nang', type=int, default=48)
    parser.add_argument('--rmin', type=float, default=1e-4)
    parser.add_argument('--length', type=float, default=20.0)
    parser.add_argument('--out', default='fem_raw.json')
    args = parser.parse_args()
    run(args.nrad, args.nang, args.rmin, args.length, args.out)
