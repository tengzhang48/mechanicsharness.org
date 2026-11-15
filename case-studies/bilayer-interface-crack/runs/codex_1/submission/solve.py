#!/usr/bin/env python3
"""Local, self-contained plane-strain calculation for PUBLIC_TASK.md.

Uses Q2 isoparametric quadrilaterals, geometrically scaled square rings
around the crack tip, and a conventional mesh away from the tip.  No
packages other than NumPy and SciPy are needed.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


NU = 1.0 / 3.0
D = np.array([[1.5, 0.75, 0.0], [0.75, 1.5, 0.0], [0.0, 0.0, 0.375]])


def exact_values():
    """Use E2=h=P(or M)=1, H=2, E1=9 throughout."""
    mu1, mu2 = 9 / (2 * (1 + NU)), 1 / (2 * (1 + NU))
    kappa = 3 - 4 * NU
    eps = np.log((kappa / mu1 + 1 / mu2) / (kappa / mu2 + 1 / mu1)) / (2 * np.pi)
    energy_factor = ((kappa + 1) / mu1 + (kappa + 1) / mu2) / (
        16 * np.cosh(np.pi * eps) ** 2
    )
    ga = (1 - NU**2) / 2 * (1 / 9 + 1 / 2 + 12 * 1.5**2 / 8)
    gb = (1 - NU**2) / 2 * (12 / 9 + 12 / 8)
    gab = (1 - NU**2) / 2 * (1 / 9 + 12 / 9 + 1 / 2 + 12 * 2.5**2 / 8)
    return {
        "eps": float(eps),
        "energy_factor": float(energy_factor),
        "G": [float(ga), float(gb), float(gab)],
        "Kabs": [float(np.sqrt(ga / energy_factor)), float(np.sqrt(gb / energy_factor))],
        "phase_difference_abs": float(np.degrees(np.arccos((gab - ga - gb) / (2 * np.sqrt(ga * gb))))),
    }


def shapes(points):
    """Q2 nodal shapes, local node order (eta=-1,0,1), then xi=-1,0,1."""
    points = np.atleast_2d(points)
    x, y = points.T
    lx = np.array([x * (x - 1) / 2, 1 - x * x, x * (x + 1) / 2]).T
    ly = np.array([y * (y - 1) / 2, 1 - y * y, y * (y + 1) / 2]).T
    dx = np.array([x - 0.5, -2 * x, x + 0.5]).T
    dy = np.array([y - 0.5, -2 * y, y + 0.5]).T
    shape = np.einsum("qi,qj->qij", ly, lx).reshape(-1, 9)
    deriv = np.stack(
        [np.einsum("qi,qj->qij", ly, dx).reshape(-1, 9),
         np.einsum("qi,qj->qij", dy, lx).reshape(-1, 9)], axis=-1
    )
    return shape, deriv


class Mesh:
    def __init__(self, nside=16, ratio=1.22, rmin=1e-7, length=24.0):
        assert nside % 2 == 0
        self.coords, self.conn, self.moduli = [], [], []
        self.lookup, self.rings, self.left_edges = {}, [], []

        # Homothetic square rings.  The innermost elements have a collapsed
        # edge at the common crack tip; Gaussian quadrature avoids the apex.
        nr = int(np.ceil(np.log(1 / rmin) / np.log(ratio)))
        radii = np.r_[0, np.geomspace(rmin, 1.0, nr + 1)]
        perimeter = []
        for j in range(nside):
            perimeter.append((1.0, -1.0 + 2 * j / nside))
        for j in range(nside):
            perimeter.append((1.0 - 2 * j / nside, 1.0))
        for j in range(nside):
            perimeter.append((-1.0, 1.0 - 2 * j / nside))
        for j in range(nside):
            perimeter.append((-1.0 + 2 * j / nside, -1.0))
        perimeter = np.array(perimeter)
        for ri, ro in zip(radii[:-1], radii[1:]):
            ring = []
            for j, a in enumerate(perimeter):
                b = perimeter[(j + 1) % len(perimeter)]
                corners = np.array([ri * a, ro * a, ro * b, ri * b])
                upper = (a[1] + b[1]) > 0
                ring.append(self.add_quad(corners, upper))
            self.rings.append((ri, ro, ring))
        self.nside = nside

        # Exterior rectangular mesh. Its central square is replaced by rings.
        step = 2 / nside
        xfar = [1.0]
        dx = step
        while xfar[-1] < length:
            xfar.append(min(length, xfar[-1] + dx))
            dx = min(dx * 1.18, 0.8)
        xs = np.r_[-np.array(xfar[:0:-1]), np.linspace(-1, 1, nside + 1), xfar[1:]]
        ys = np.linspace(-2, 1, 3 * nside // 2 + 1)
        for xa, xb in zip(xs[:-1], xs[1:]):
            for ya, yb in zip(ys[:-1], ys[1:]):
                if xa >= -1 and xb <= 1 and ya >= -1:
                    continue
                corners = np.array([[xa, ya], [xb, ya], [xb, yb], [xa, yb]])
                ie = self.add_quad(corners, ya >= 0)
                if xa == -length:
                    self.left_edges.append(ie)
        self.coords = np.array(self.coords)
        self.conn = np.array(self.conn, dtype=np.int32)
        self.moduli = np.array(self.moduli)
        self.edofs = (2 * self.conn[:, :, None] + np.arange(2)).reshape(-1, 18)
        self.length = length

    @staticmethod
    def key(x, y, upper):
        # Rounding only serves to identify the same geometric node built
        # through two adjacent elements. Stored coordinates are unrounded.
        x, y = round(float(x), 14), round(float(y), 14)
        crack_face = int(upper) if y == 0.0 and x < 0 else 0
        return x, y, crack_face

    def node(self, x, y, upper):
        key = self.key(x, y, upper)
        if key not in self.lookup:
            self.lookup[key] = len(self.coords)
            self.coords.append((float(x), float(y)))
        return self.lookup[key]

    def add_quad(self, corners, upper):
        indices = []
        a, b, c, d = corners
        for eta in (0, 0.5, 1):
            for xi in (0, 0.5, 1):
                pt = (1 - eta) * ((1 - xi) * a + xi * b) + eta * ((1 - xi) * d + xi * c)
                indices.append(self.node(*pt, upper))
        ie = len(self.conn)
        self.conn.append(indices)
        self.moduli.append(9.0 if upper else 1.0)
        return ie


def geometry(xy, deriv):
    jac = np.einsum("eni,qnj->eqij", xy, deriv)
    det = np.linalg.det(jac)
    if np.any(det <= 0):
        raise ValueError("Nonpositive element Jacobian")
    grad = np.einsum("qnj,eqji->eqni", deriv, np.linalg.inv(jac))
    return grad, det, jac


def assemble(mesh, nq=4):
    ndof = 2 * len(mesh.coords)
    ne = len(mesh.conn)
    xg, wg = leggauss(nq)
    points = np.array([(x, y) for y in xg for x in xg])
    weights = np.outer(wg, wg).ravel()
    _, deriv = shapes(points)
    # Retain all element entries until converting COO to CSR.
    values = np.empty((ne, 18, 18))
    for start in range(0, ne, 512):
        sl = slice(start, min(start + 512, ne))
        grad, det, _ = geometry(mesh.coords[mesh.conn[sl]], deriv)
        nx, ny = grad[..., 0], grad[..., 1]
        weights_e = det * weights * mesh.moduli[sl, None]
        xx = np.einsum("eqi,eqj,eq->eij", nx, nx, weights_e, optimize=True)
        yy = np.einsum("eqi,eqj,eq->eij", ny, ny, weights_e, optimize=True)
        xy = np.einsum("eqi,eqj,eq->eij", nx, ny, weights_e, optimize=True)
        values[sl, 0::2, 0::2] = 1.5 * xx + 0.375 * yy
        values[sl, 1::2, 1::2] = 0.375 * xx + 1.5 * yy
        values[sl, 0::2, 1::2] = 0.75 * xy + 0.375 * xy.transpose(0, 2, 1)
        values[sl, 1::2, 0::2] = values[sl, 0::2, 1::2].transpose(0, 2, 1)
    rows = np.broadcast_to(mesh.edofs[:, :, None], values.shape).ravel()
    cols = np.broadcast_to(mesh.edofs[:, None, :], values.shape).ravel()
    stiff = coo_matrix((values.ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()
    force = np.zeros((ndof, 2))
    # On the left end the outward traction is (-sigma11,0).
    sg, sd = shapes(np.column_stack([-np.ones(nq), xg]))
    for ie in mesh.left_edges:
        xy = mesh.coords[mesh.conn[ie]]
        pts = sg @ xy
        jac_edge = np.linalg.norm(sd[:, :, 1] @ xy, axis=1)
        y = pts[:, 1]
        if mesh.moduli[ie] == 9:
            stresses = np.column_stack([-np.ones(nq), -12 * (y - 0.5)])
        else:
            stresses = np.column_stack([0.5 + 2.25 * (y + 1), 1.5 * (y + 1)])
        load = sg.T @ (-stresses * (wg * jac_edge)[:, None])
        np.add.at(force, mesh.edofs[ie, 0::2], load)
    return stiff, force


def field(mesh, displacement, ie, points):
    """Returns displacements gradients and stresses, each with both cases."""
    _, deriv = shapes(points)
    grad, _, jac = geometry(mesh.coords[mesh.conn[ie]][None], deriv)
    ue = displacement[mesh.edofs[ie]].reshape(9, 2, 2)
    # q, displacement component, spatial derivative, load case
    du = np.einsum("qnj,nic->qijc", grad[0], ue)
    strain = np.stack([du[:, 0, 0], du[:, 1, 1], du[:, 0, 1] + du[:, 1, 0]], axis=1)
    stress = mesh.moduli[ie] * np.einsum("ij,qjc->qic", D, strain)
    return du, strain, stress, jac[0]


def extract(mesh, displacement):
    exact = exact_values()
    eps = exact["eps"]
    # Crack opening: [u2]+i[u1] = 4 S / cosh(pi eps) K sqrt(r/(2pi))
    #                           * r^(i eps)/(1+2i eps), S=sum(1/Ek').
    compliance_sum = (1 - NU**2) * (1 + 1 / 9)
    cod_factor = np.sqrt(2 * np.pi) * (1 + 2j * eps) * np.cosh(np.pi * eps) / (4 * compliance_sum)
    cod_rows, traction_rows, j_rows = [], [], []
    for ri, ro, ring in mesh.rings[1:]:
        # Use the ring outer nodes, which belong to both adjoining elements.
        nt = mesh.lookup[mesh.key(-ro, 0.0, True)]
        nb = mesh.lookup[mesh.key(-ro, 0.0, False)]
        delta = displacement[2 * nt + 1] - displacement[2 * nb + 1] + 1j * (
            displacement[2 * nt] - displacement[2 * nb]
        )
        kc = cod_factor * delta / (ro ** (0.5 + 1j * eps))
        cod_rows.append([ro, *kc])
        # Tractions on the interface ahead, evaluated from either material.
        for xi in (-1 / np.sqrt(3), 0.0, 1 / np.sqrt(3)):
            r = (ri + ro) / 2 + (ro - ri) * xi / 2
            ku = []
            for ie, eta in [(ring[mesh.nside // 2], -1.0), (ring[mesh.nside // 2 - 1], 1.0)]:
                _, _, stress, _ = field(mesh, displacement, ie, [(xi, eta)])
                k = (stress[0, 1] + 1j * stress[0, 2]) * np.sqrt(2 * np.pi * r) * r ** (-1j * eps)
                ku.extend(k)
            traction_rows.append([r, xi, *ku])

        # J integral on the middle square of every ring.
        if 1e-4 < ro < 0.9:
            zg, zw = leggauss(5)
            points = np.column_stack([np.zeros(5), zg])
            jvalue = np.zeros(2)
            for ie in ring:
                du, strain, stress, jac = field(mesh, displacement, ie, points)
                dx, dy = jac[:, 0, 1], jac[:, 1, 1]
                energy = 0.5 * np.sum(strain * stress, axis=1)
                integral = energy * dy[:, None] - (
                    stress[:, 0] * dy[:, None] - stress[:, 2] * dx[:, None]
                ) * du[:, 0, 0] - (
                    stress[:, 2] * dy[:, None] - stress[:, 1] * dx[:, None]
                ) * du[:, 1, 0]
                jvalue += zw @ integral
            j_rows.append([(ri + ro) / 2, *jvalue])
    return np.array(cod_rows), np.array(traction_rows), np.array(j_rows)


def complex_list(z):
    return [[float(v.real), float(v.imag)] for v in z]


def estimate(cod, lower=1e-4, upper=1e-2):
    sel = (cod[:, 0].real >= lower) & (cod[:, 0].real <= upper)
    r = cod[sel, 0].real
    # Leading Williams field plus the next, O(r), correction to K.
    design = np.column_stack([np.ones(len(r)), r, r**2])
    return np.linalg.lstsq(design, cod[sel, 1:], rcond=None)[0][0]


def run(args):
    start = time.time()
    mesh = Mesh(args.nside, args.ratio, args.rmin, args.length)
    print(f"mesh: {len(mesh.coords)} nodes, {len(mesh.conn)} elements; {time.time()-start:.2f}s", flush=True)
    stiff, force = assemble(mesh, args.quadrature)
    n = len(mesh.coords)
    loads = force.reshape(n, 2, 2)
    moment = np.einsum("n,nc->c", mesh.coords[:, 0], loads[:, 1]) - np.einsum(
        "n,nc->c", mesh.coords[:, 1], loads[:, 0]
    )
    print(f"force balance: {np.sum(loads, axis=0)}; moment: {moment}", flush=True)
    p = mesh.lookup[mesh.key(args.length, 0, True)]
    q = mesh.lookup[mesh.key(args.length, 1, True)]
    fixed = [2 * p, 2 * p + 1, 2 * q]
    free = np.setdiff1d(np.arange(2 * n), fixed)
    print(f"assembled: {time.time()-start:.2f}s, nnz={stiff.nnz}", flush=True)
    factor = splu(stiff[free][:, free].tocsc())
    displacement = np.zeros_like(force)
    displacement[free] = factor.solve(force[free])
    residual = stiff @ displacement - force
    print(f"solved: {time.time()-start:.2f}s; max residual {np.max(np.abs(residual)):.3g}", flush=True)
    cod, tractions, jvals = extract(mesh, displacement)
    k = estimate(cod)
    phases = np.degrees(np.angle(np.r_[k, k.sum()]))
    summary = {
        "parameters": vars(args),
        "nodes": n,
        "elements": len(mesh.conn),
        "exact": exact_values(),
        "K_cod": complex_list(k),
        "Kabs_cod": np.abs(k).tolist(),
        "phases_deg_cod": phases.tolist(),
        "residual_max": float(np.max(np.abs(residual))),
        "force_balance": np.sum(loads, axis=0).tolist(),
        "moment_balance": moment.tolist(),
        "elapsed_seconds": time.time() - start,
    }
    for lo, hi in [(1e-5, 1e-3), (1e-4, 1e-2), (1e-3, 0.03)]:
        kk = estimate(cod, lo, hi)
        summary[f"fit_{lo:g}_{hi:g}"] = {"K": complex_list(kk), "phases": np.degrees(np.angle(np.r_[kk, kk.sum()])).tolist()}
    prefix = Path(args.output)
    prefix.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    np.savez_compressed(prefix.with_suffix(".npz"), cod=cod, tractions=tractions, j=jvals)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nside", type=int, default=16)
    parser.add_argument("--ratio", type=float, default=1.22)
    parser.add_argument("--rmin", type=float, default=1e-7)
    parser.add_argument("--length", type=float, default=24.0)
    parser.add_argument("--quadrature", type=int, default=4)
    parser.add_argument("--output", default="fem_default")
    run(parser.parse_args())
