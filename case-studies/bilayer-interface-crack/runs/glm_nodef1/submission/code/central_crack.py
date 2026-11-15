"""Central interface crack (length 2a, tips at 0 and -2a) in a long strip
under remote uniform tension s11 = sigma_inf at both ends.  Compare K at the
right tip with the classical Rice-Sih form and with G-based predictions.
"""
import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import (plane_strain_D, assemble, solve, geom_breaks, q8_shape,
                 elem_stress_at, end_loads)
from fits import fit_series

nu = 1/3
E1, E2 = 9.0, 1.0
D = [plane_strain_D(E1, nu), plane_strain_D(E2, nu)]
h, H = 1.0, 2.0
a = 1.0
L = 40.0
eps = np.log(2/3)/(2*np.pi)
S = (1-nu**2)/E1 + (1-nu**2)/E2
Eb1, Eb2 = E1/(1-nu**2), E2/(1-nu**2)
Sb = 1/Eb1 + 1/Eb2


def build_X(d1):
    A = geom_breaks(a, d1, 1.4, 0.4)          # 0..a fine at 0
    C = -A[::-1]                               # -a..0 fine at 0
    B = -2*a + A                                # -2a..-a fine at -2a
    Dx = -2*a - geom_breaks(L-2*a, d1, 1.6, 1.5)  # -L..-2a fine at -2a
    E = geom_breaks(L, d1, 1.6, 1.5)            # 0..L fine at 0
    X = np.unique(np.concatenate([Dx, B, C, E]))
    return X


class CentralCrack(fem.BilayerCrack):
    def __init__(self, h, H, L, d1, d2, a):
        # patch: build custom X then call a modified constructor body
        self.a = a
        Yu = geom_breaks(h, d2, 1.5, 0.35)
        Yl = -geom_breaks(H, d2, 1.5, 0.35)[::-1][:-1]
        self.Y = np.concatenate([Yl, Yu])
        self.j0 = len(Yl)
        X = build_X(d1)
        self.X = X
        self.it = int(np.where(X == 0)[0][0])
        nx, ny = len(X), len(self.Y)
        nid = np.arange(nx*ny).reshape(nx, ny)
        hbase = nx*ny
        hnid = hbase + np.arange((nx-1)*ny).reshape(nx-1, ny)
        vbase = hbase + (nx-1)*ny
        vnid = vbase + np.arange(nx*(ny-1)).reshape(nx, ny-1)
        nn = vbase + nx*(ny-1)
        coords = np.zeros((nn, 2))
        for i in range(nx):
            for j in range(ny):
                coords[nid[i, j]] = (X[i], self.Y[j])
        for i in range(nx-1):
            for j in range(ny):
                coords[hnid[i, j]] = (0.5*(X[i]+X[i+1]), self.Y[j])
        for i in range(nx):
            for j in range(ny-1):
                coords[vnid[i, j]] = (X[i], 0.5*(self.Y[j]+self.Y[j+1]))
        remap = np.arange(nn)
        dup = []
        new = nn
        for i in range(nx):
            if -2*a < X[i] < 0:
                b = nid[i, self.j0]
                dup.append((b, new, X[i])); remap[b] = new; new += 1
        for i in range(nx-1):
            xm = 0.5*(X[i]+X[i+1])
            if -2*a < xm < 0:
                b = hnid[i, self.j0]
                dup.append((b, new, xm)); remap[b] = new; new += 1
        if dup:
            coords = np.vstack([coords, [(x, 0.0) for (_, _, x) in dup]])
        self.coords = coords
        self.nid, self.hnid, self.vnid = nid, hnid, vnid
        self.nx, self.ny = nx, ny
        elems, mats = [], []
        for j in range(ny-1):
            upper = j >= self.j0
            for i in range(nx-1):
                e = [nid[i, j], nid[i+1, j], nid[i+1, j+1], nid[i, j+1],
                     hnid[i, j], vnid[i+1, j], hnid[i, j+1], vnid[i, j]]
                if not upper:
                    e = [remap[k] for k in e]
                elems.append(e)
                mats.append(1 if upper else 2)
        self.elems = np.array(elems)
        self.mats = np.array(mats)
        self.elem_of_cell = self.elems.reshape(ny-1, nx-1, 8)
        self.pairs = [(b, d, x) for (b, d, x) in dup]
        edges = []
        for j in range(ny-1):
            n = [nid[0, j], nid[0, j+1], vnid[0, j]]
            if j < self.j0:
                n = [remap[k] for k in n]
            edges.append((n, self.Y[j], self.Y[j+1], 1 if j >= self.j0 else 2))
        self.left_edges = edges
        redges = []
        for j in range(ny-1):
            n = [nid[nx-1, j], nid[nx-1, j+1], vnid[nx-1, j]]
            redges.append((n, self.Y[j], self.Y[j+1], 1 if j >= self.j0 else 2))
        self.right_edges = redges


def end_loads_uniform(mesh, sigma, both=False):
    f = np.zeros(2*len(mesh.coords))
    pairs = ((mesh.left_edges, -1.0), (mesh.right_edges, +1.0)) if both \
        else ((mesh.left_edges, -1.0),)
    for edges, sgn in pairs:
        for (n, yb, yt, mat) in edges:
            half = 0.5*(yt-yb)
            for (t, w) in ((-np.sqrt(3/5), 5/9), (0.0, 8/9), (np.sqrt(3/5), 5/9)):
                Na = 0.5*t*(t-1); Nb = 0.5*t*(t+1); Nm = 1-t*t
                fx = sgn*sigma*half*w
                f[2*n[0]] += Na*fx
                f[2*n[1]] += Nb*fx
                f[2*n[2]] += Nm*fx
    return f


mesh = CentralCrack(h, H, L, 3e-4, 2e-4, a)
f = end_loads_uniform(mesh, 1.0)
u = solve(mesh, D, f)

rs = np.geomspace(2e-3, 5e-2, 40)
T = fem.interface_tractions(mesh, D, u, rs)
for (rmin, rmax) in [(3e-3, 5e-2), (1e-3, 2e-2)]:
    m = (rs >= rmin) & (rs <= rmax)
    y = T[m]*np.sqrt(2*np.pi*rs[m])*np.exp(-1j*eps*np.log(rs[m]))
    K = fit_series(rs[m], y, (0.0, 0.5, 1.0))[0]
    print(f"window {rmin:.0e}-{rmax:.0e}: K = {K:.6f}, |K|/sqrt(pi a) = {abs(K)/np.sqrt(np.pi):.5f}, "
          f"phase = {np.degrees(np.angle(K)):+.3f} deg")

print("\ncandidates (sigma=1, a=1):")
print(f"  K = (1+2i eps) cosh(pi eps) sqrt(pi a): "
      f"{(1+2j*eps)*np.cosh(np.pi*eps)*np.sqrt(np.pi):.5f} "
      f"phase {np.degrees(np.angle((1+2j*eps)*np.cosh(np.pi*eps)*np.sqrt(np.pi))):+.3f}")
print(f"  K = conj(1+2i eps) cosh(pi eps) sqrt(pi a): "
      f"{np.conj(1+2j*eps)*np.cosh(np.pi*eps)*np.sqrt(np.pi):.5f} "
      f"phase {np.degrees(np.angle(np.conj(1+2j*eps))):+.3f}")
print(f"  K = cosh(pi eps) sqrt(pi a): {np.cosh(np.pi*eps)*np.sqrt(np.pi):.5f}")
