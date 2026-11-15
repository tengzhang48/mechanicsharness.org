"""Central-crack benchmark, clean mesh: fine at tips 0 and -2a, coarse between
and beyond."""
import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import (plane_strain_D, solve, geom_breaks, interface_tractions,
                 crack_opening, BilayerCrack)
from fits import fit_series

nu = 1/3
E1, E2 = 9.0, 1.0
D = [plane_strain_D(E1, nu), plane_strain_D(E2, nu)]
h, H, a, L = 1.0, 2.0, 1.0, 40.0
eps = np.log(2/3)/(2*np.pi)

# X construction: fine at 0 and at -2a
S1 = geom_breaks(a, 4e-4, 1.4, 0.4)             # 0..a, fine at 0
S2 = geom_breaks(L-a, 0.4, 1.5, 2.0) + a         # a..L, coarse
pos = np.concatenate([S1, S2[1:]])
neg1 = -S1[::-1]                                  # -a..0 fine at 0
neg2 = (-2*a + np.concatenate([S1, S2[1:]]))[::-1]  # -L..-2a? build -2a..-a then... 
# simpler: mirror: right side of left tip:
Xr = pos                                          # 0..L
Xl = -pos                                         # -L..0 (fine at 0)
Xl2 = -(2*a - S1[::-1])                           # = -2a + S1: -2a..-a fine at -2a
# assemble: [-L..-2a] coarse-fine-at(-2a) + [-2a..-a] fine at -2a + [-a..0] fine at 0 + [0..L]
segA = (-2*a - geom_breaks(L-a, 0.4, 1.5, 2.0))[::-1]   # -L..-2a coarse, fine at -2a
X = np.unique(np.concatenate([segA, Xl2, neg1[1:], Xr]))
print("nx:", len(X), "min dx in (-2a,-a) region:",
      np.diff(X[(X > -2*a) & (X < 0)]).min())


class CC(BilayerCrack):
    pass


# reuse BilayerCrack-style construction with custom X and crack range:
class Central2(BilayerCrack):
    def __init__(self, h, H, Lx, d1, d2, a, X):
        Yu = geom_breaks(h, d2, 1.5, 0.35)
        Yl = -geom_breaks(H, d2, 1.5, 0.35)[::-1][:-1]
        self.Y = np.concatenate([Yl, Yu])
        self.j0 = len(Yl)
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
        coords[nid.ravel(), 0] = np.repeat(X, ny)
        coords[nid.ravel(), 1] = np.tile(self.Y, nx)
        hm_x = 0.5*(X[:-1]+X[1:])
        coords[hnid.ravel(), 0] = np.repeat(hm_x, ny)
        coords[hnid.ravel(), 1] = np.tile(self.Y, nx-1)
        coords[vnid.ravel(), 0] = np.repeat(X, ny-1)
        coords[vnid.ravel(), 1] = np.tile(0.5*(self.Y[:-1]+self.Y[1:]), nx)
        remap = np.arange(nn)
        dup = []
        new = nn
        for i in range(nx):
            if -2*a < X[i] < 0:
                b = int(nid[i, self.j0]); dup.append((b, new, X[i])); remap[b] = new; new += 1
        for i in range(nx-1):
            xm = 0.5*(X[i]+X[i+1])
            if -2*a < xm < 0:
                b = int(hnid[i, self.j0]); dup.append((b, new, xm)); remap[b] = new; new += 1
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
        self.pairs = dup
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


mesh = Central2(h, H, L, 4e-4, 2e-4, a, X)
f = np.zeros(2*len(mesh.coords))
for edges, sgn in ((mesh.left_edges, -1.0), (mesh.right_edges, +1.0)):
    for (n, yb, yt, mat) in edges:
        half = 0.5*(yt-yb)
        for (t, w) in ((-np.sqrt(3/5), 5/9), (0.0, 8/9), (np.sqrt(3/5), 5/9)):
            Na = 0.5*t*(t-1); Nb = 0.5*t*(t+1); Nm = 1-t*t
            fx = sgn*1.0*half*w
            f[2*n[0]] += Na*fx
            f[2*n[1]] += Nb*fx
            f[2*n[2]] += Nm*fx
u = solve(mesh, D, f)
rd, dd = crack_opening(mesh, u)
m = np.argmin(np.abs(rd-1.0))
print(f"max|u|={np.abs(u).max():.4f}  delta(x=-1)={dd[m]:.5f}")
rs = np.geomspace(3e-3, 5e-2, 30)
T = interface_tractions(mesh, D, u, rs)
for (rmin, rmax) in [(3e-3, 5e-2), (1e-3, 2e-2)]:
    mm = (rs >= rmin) & (rs <= rmax)
    y = T[mm]*np.sqrt(2*np.pi*rs[mm])*np.exp(-1j*eps*np.log(rs[mm]))
    K = fit_series(rs[mm], y, (0.0, 0.5, 1.0))[0]
    print(f"win {rmin:.0e}-{rmax:.0e}: K={K:.5f} |K|/sqrt(pi a)={abs(K)/np.sqrt(np.pi):.5f} "
          f"phase={np.degrees(np.angle(K)):+.3f}")
print("candidate (1+2 i eps)cosh(pi eps)sqrt(pi):",
      (1+2j*eps)*np.cosh(np.pi*eps)*np.sqrt(np.pi))
