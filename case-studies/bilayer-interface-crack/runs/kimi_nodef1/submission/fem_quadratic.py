#!/usr/bin/env python3
"""Quadratic (Q9/T6) version of the bilayer crack finite-element model."""
from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import fem_solve as base


@dataclass
class QuadraticMesh:
    xy: np.ndarray
    top_elements: np.ndarray
    bottom_elements: np.ndarray
    top_triangles: np.ndarray
    bottom_triangles: np.ndarray
    top_outer: np.ndarray
    bottom_outer: np.ndarray
    edge_mid: dict
    crack_top: np.ndarray
    crack_bottom: np.ndarray
    crack_r: np.ndarray
    right_tip: int
    right_top_corner: int


def elevate(m: base.Mesh) -> QuadraticMesh:
    coordinates = m.xy.tolist()
    edge_mid = {}

    def midpoint(a, b):
        key = (min(int(a), int(b)), max(int(a), int(b)))
        if key not in edge_mid:
            edge_mid[key] = len(coordinates)
            coordinates.append(((m.xy[a] + m.xy[b]) / 2.0).tolist())
        return edge_mid[key]

    top_q, bottom_q = [], []
    for source, destination in [(m.top_elements, top_q),
                                (m.bottom_elements, bottom_q)]:
        for corners in source:
            a, b, c, d = map(int, corners)
            sides = [midpoint(a, b), midpoint(b, c),
                     midpoint(c, d), midpoint(d, a)]
            center = len(coordinates)
            coordinates.append(np.mean(m.xy[corners], axis=0).tolist())
            destination.append([a, b, c, d] + sides + [center])
    top_t, bottom_t = [], []
    for source, destination in [(m.top_triangles, top_t),
                                (m.bottom_triangles, bottom_t)]:
        for corners in source:
            a, b, c = map(int, corners)
            destination.append([a, b, c, midpoint(a, b), midpoint(b, c),
                                midpoint(c, a)])

    ct, cb, cr = [], [], []
    for j in range(len(m.crack_r) - 1):
        ta, tb = int(m.crack_top[j]), int(m.crack_top[j + 1])
        ba, bb = int(m.crack_bottom[j]), int(m.crack_bottom[j + 1])
        if j == 0:
            ct.append(ta)
            cb.append(ba)
            cr.append(m.crack_r[j])
        ct.extend([midpoint(ta, tb), tb])
        cb.extend([midpoint(ba, bb), bb])
        tm, bm = edge_mid[(min(ta, tb), max(ta, tb))], \
                 edge_mid[(min(ba, bb), max(ba, bb))]
        cr.extend([abs(coordinates[tm][0]), m.crack_r[j + 1]])
        assert abs(coordinates[tm][0] - coordinates[bm][0]) < 1e-12
    return QuadraticMesh(np.asarray(coordinates), np.asarray(top_q),
                         np.asarray(bottom_q), np.asarray(top_t),
                         np.asarray(bottom_t), m.top_outer, m.bottom_outer,
                         edge_mid, np.asarray(ct), np.asarray(cb),
                         np.asarray(cr), m.right_tip, m.right_top_corner)


def q9_stiffness(xy, elements, D):
    gp = math.sqrt(3.0 / 5.0)
    gauss = [-gp, 0.0, gp]
    weights = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
    # Corners, side midpoints, center in the conventional Q9 ordering.
    locations = [(-1, -1), (1, -1), (1, 1), (-1, 1),
                 (0, -1), (1, 0), (0, 1), (-1, 0), (0, 0)]
    pxy = xy[elements]
    ke = np.zeros((len(elements), 18, 18))
    for xi, wx in zip(gauss, weights):
        Lx = np.array([xi * (xi - 1) / 2, 1 - xi * xi,
                       xi * (xi + 1) / 2])
        dx = np.array([xi - 0.5, -2 * xi, xi + 0.5])
        for eta, wy in zip(gauss, weights):
            Ly = np.array([eta * (eta - 1) / 2, 1 - eta * eta,
                           eta * (eta + 1) / 2])
            dy = np.array([eta - 0.5, -2 * eta, eta + 0.5])
            dN = np.zeros((9, 2))
            for g, (a, b) in enumerate(locations):
                ia, ib = a + 1, b + 1
                dN[g] = [dx[ia] * Ly[ib], Lx[ia] * dy[ib]]
            jac = np.einsum('gi,egj->eij', dN, pxy)
            det = np.linalg.det(jac)
            if np.any(det <= 0):
                raise ValueError(f'Nonpositive Q9 Jacobian {det.min()}')
            dxy = np.einsum('eij,gj->egi', np.linalg.inv(jac), dN)
            B = np.zeros((len(elements), 3, 18))
            B[:, 0, 0::2] = dxy[:, :, 0]
            B[:, 1, 1::2] = dxy[:, :, 1]
            B[:, 2, 0::2] = dxy[:, :, 1]
            B[:, 2, 1::2] = dxy[:, :, 0]
            ke += wx * wy * np.einsum('eki,kl,elj,e->eij', B, D, B, det)
    return ke


def t6_stiffness(xy, elements, D):
    p = xy[elements]
    x1, y1 = p[:, 0, 0], p[:, 0, 1]
    x2, y2 = p[:, 1, 0], p[:, 1, 1]
    x3, y3 = p[:, 2, 0], p[:, 2, 1]
    det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if np.any(det <= 0):
        raise ValueError('Nonpositive T6 Jacobian')
    area = det / 2
    # Constant metric derivatives.
    grads = np.zeros((len(elements), 3, 2))
    grads[:, 0] = np.stack((y2 - y3, x3 - x2), axis=1) / det[:, None]
    grads[:, 1] = np.stack((y3 - y1, x1 - x3), axis=1) / det[:, None]
    grads[:, 2] = np.stack((y1 - y2, x2 - x1), axis=1) / det[:, None]
    ke = np.zeros((len(elements), 12, 12))
    for a, b in [(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)]:
        c = 1 - a - b
        # Node order is corner 1,2,3 then midpoint 12,23,31. The
        # following are partial derivatives with respect to L1, L2, L3.
        darea = np.array([[4*a-1, 0, 0, 4*b, 0, 4*c],
                          [0, 4*b-1, 0, 4*a, 4*c, 0],
                          [0, 0, 4*c-1, 0, 4*b, 4*a]])
        dNxy = np.einsum('an,eap->enp', darea, grads)
        B = np.zeros((len(elements), 3, 12))
        B[:, 0, 0::2] = dNxy[:, :, 0]
        B[:, 1, 1::2] = dNxy[:, :, 1]
        B[:, 2, 0::2] = dNxy[:, :, 1]
        B[:, 2, 1::2] = dNxy[:, :, 0]
        ke += (area / 3)[:, None, None] * np.einsum('eki,kl,elj->eij', B, D, B)
    return ke


def assemble(mesh, E1=9.0, E2=1.0, nu=1/3):
    n = len(mesh.xy)
    blocks = []
    for els, E in [(mesh.top_elements, E1), (mesh.bottom_elements, E2)]:
        blocks.append((els, q9_stiffness(mesh.xy, els,
                                         base.plane_strain_D(E, nu))))
    for els, E in [(mesh.top_triangles, E1),
                   (mesh.bottom_triangles, E2)]:
        blocks.append((els, t6_stiffness(mesh.xy, els,
                                         base.plane_strain_D(E, nu))))
    rows, cols, vals = [], [], []
    for els, ke in blocks:
        edof = np.empty((len(els), els.shape[1] * 2), dtype=np.int64)
        edof[:, 0::2] = 2 * els
        edof[:, 1::2] = 2 * els + 1
        rows.append(np.repeat(edof, edof.shape[1], axis=1).ravel())
        cols.append(np.tile(edof, (1, edof.shape[1])).ravel())
        vals.append(ke.ravel())
    return sp.coo_matrix((np.concatenate(vals),
                          (np.concatenate(rows), np.concatenate(cols))),
                         shape=(2*n, 2*n)).tocsr()


def consistent_loads(mesh, length, h=1.0, H=2.0):
    n = len(mesh.xy)
    F = np.zeros((2*n, 2))
    gp = math.sqrt(3/5)
    gauss = [(-gp, 5/9), (0, 8/9), (gp, 5/9)]
    # Lagrange basis at endpoints and midpoint of a quadratic edge.
    for outer, upper in [(mesh.top_outer, True),
                         (mesh.bottom_outer, False)]:
        for a, b in zip(outer[:-1], outer[1:]):
            a, b = int(a), int(b)
            xa, ya = mesh.xy[a]
            xb, yb = mesh.xy[b]
            if abs(xa + length) > 1e-9 or abs(xb + length) > 1e-9 or abs(xa-xb) > 1e-9:
                continue
            mid = mesh.edge_mid[(min(a,b), max(a,b))]
            for P, M, column in [(1.0, 0.0, 0), (0.0, 1.0, 1)]:
                for xi, w in gauss:
                    s = (xi+1)/2
                    y = ya*(1-s) + yb*s
                    if upper:
                        stress = -P/h - 12*M*(y-h/2)/h**3
                    else:
                        M2 = M + P*(h+H)/2
                        stress = P/H + 12*M2*(y+H/2)/H**3
                    tx = -stress
                    jac = abs(yb-ya)/2
                    N = np.array([xi*(xi-1)/2, xi*(xi+1)/2, 1-xi*xi])
                    # a and b correspond to natural -1 and +1 respectively.
                    for node, Ni in zip([a,b,mid], N):
                        F[2*node,column] += Ni*tx*w*jac
    return F


def run(nrad=100, nang=24, rmin=1e-4, length=20, out='quadratic_raw.json'):
    start = time.time()
    first = base.build_mesh(nrad, nang, rmin, length)
    mesh = elevate(first)
    print(f'Q9 mesh: {len(mesh.xy)} nodes, {len(mesh.top_elements)+len(mesh.bottom_elements)} quads')
    K = assemble(mesh)
    print(f'stiffness: {K.nnz} entries at {time.time()-start:.1f}s')
    F = consistent_loads(mesh, length)
    fixed = np.array([2*mesh.right_tip, 2*mesh.right_tip+1,
                      2*mesh.right_top_corner])
    free = np.setdiff1d(np.arange(2*len(mesh.xy)), fixed)
    U = np.zeros_like(F)
    Kfree = K[free][:,free].tocsr()
    # Geometric grading gives diagonal entries spanning many orders of
    # magnitude. Symmetric Jacobi equilibration improves both sparse pivoting
    # and the accuracy of the tiny crack-face displacements.
    scale = 1.0 / np.sqrt(Kfree.diagonal())
    Kscaled = Kfree.multiply(scale[:, None]).multiply(scale[None, :]).tocsr()
    U[free] = scale[:, None] * spla.spsolve(Kscaled,
                                             F[free] * scale[:, None])
    residual = K @ U - F
    relative = np.linalg.norm(residual[free], axis=0) / np.linalg.norm(F[free], axis=0)
    print(f'solved at {time.time()-start:.1f}s; free residuals {relative}')
    nu = 1/3
    kappa = 3-4*nu
    mu1, mu2 = 9/(2*(1+nu)), 1/(2*(1+nu))
    beta = ((mu1-mu2)*(kappa-1))/((mu1+mu2)*(kappa+1))
    eps = math.log((1-beta)/(1+beta))/(2*math.pi)
    D = (kappa+1)/mu1 + (kappa+1)/mu2
    delta = []
    for col in range(2):
        ut = U[:,col].reshape(-1,2)[mesh.crack_top]
        ub = U[:,col].reshape(-1,2)[mesh.crack_bottom]
        d = ut-ub
        delta.append(d[:,1]+1j*d[:,0])
    result = {'parameters': {'nrad':nrad,'nang':nang,'rmin':rmin,
                             'length':length,'nodes':len(mesh.xy),
                             'order':2,'beta':beta,'epsilon':eps,'D':D},
              'r':mesh.crack_r.tolist()}
    for name, d in zip(['A','B'],delta):
        result['delta'+name+'_real'] = d.real.tolist()
        result['delta'+name+'_imag'] = d.imag.tolist()
    with open(out,'w') as f:
        json.dump(result,f)
    print(f'wrote {out}; elapsed {time.time()-start:.1f}s')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--nrad',type=int,default=100)
    p.add_argument('--nang',type=int,default=24)
    p.add_argument('--rmin',type=float,default=1e-4)
    p.add_argument('--length',type=float,default=20)
    p.add_argument('--out',default='quadratic_raw.json')
    a=p.parse_args()
    run(a.nrad,a.nang,a.rmin,a.length,a.out)
