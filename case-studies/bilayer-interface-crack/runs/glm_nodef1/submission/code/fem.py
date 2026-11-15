"""Q8 (8-node quad) plane-strain FEM for the bilayer interface-crack problem.

Geometry: layer 1 (material 1) occupies 0 < x2 < h, layer 2 (material 2)
occupies -H < x2 < 0.  Interface crack along x1 < 0 with tip at origin.
Outer surfaces x2 = h, x2 = -H traction free.  End x1 = +L traction free.
End x1 = -L loaded with the Saint-Venant stress distribution given in the
task (consistent nodal forces).  Loads are self-equilibrated so the far
ahead region is stress free.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

# ---------------------------------------------------------------- Q8 element
# node order: corners (0,1,2,3) CCW from (-1,-1), then midsides 4..7 on
# edges 0-1, 1-2, 2-3, 3-0.
_GL = np.sqrt(3.0 / 5.0)
_GP = []
for a in (-_GL, 0.0, _GL):
    for b in (-_GL, 0.0, _GL):
        w = (5/9 if a != 0 else 8/9) * (5/9 if b != 0 else 8/9)
        _GP.append((a, b, w))


def q8_shape(xi, eta):
    N = np.array([
        0.25*(1-xi)*(1-eta)*(-xi-eta-1),
        0.25*(1+xi)*(1-eta)*( xi-eta-1),
        0.25*(1+xi)*(1+eta)*( xi+eta-1),
        0.25*(1-xi)*(1+eta)*(-xi+eta-1),
        0.5*(1-xi**2)*(1-eta),
        0.5*(1+xi)*(1-eta**2),
        0.5*(1-xi**2)*(1+eta),
        0.5*(1-xi)*(1-eta**2)])
    dNxi = np.array([
        0.25*(1-eta)*(2*xi+eta),
        0.25*(1-eta)*(2*xi-eta),
        0.25*(1+eta)*(2*xi+eta),
        0.25*(1+eta)*(2*xi-eta),
        -xi*(1-eta),
        0.5*(1-eta**2),
        -xi*(1+eta),
        -0.5*(1-eta**2)])
    dNeta = np.array([
        0.25*(1-xi)*(2*eta+xi),
        0.25*(1+xi)*(2*eta-xi),
        0.25*(1+xi)*(xi+2*eta),
        0.25*(1-xi)*(-xi+2*eta),
        -0.5*(1-xi**2),
        -eta*(1+xi),
        0.5*(1-xi**2),
        -eta*(1-xi)])
    return N, dNxi, dNeta


def plane_strain_D(E, nu):
    lam = E*nu/((1+nu)*(1-2*nu))
    mu = E/(2*(1+nu))
    return np.array([[lam+2*mu, lam, 0.0],
                     [lam, lam+2*mu, 0.0],
                     [0.0, 0.0, mu]])


def assemble(coords, elems, mats, D_list):
    """Sparse stiffness (2 dof/node). elems (n,8), mats (n,) int indices."""
    n = len(elems)
    xy = coords[elems]                      # (n,8,2)
    Ke = np.zeros((n, 16, 16))
    Ds = np.array(D_list)                   # (nm,3,3); mats are 1-based
    D_e = Ds[mats-1]                        # (n,3,3)
    for (xi, eta, w) in _GP:
        N, dxi, deta = q8_shape(xi, eta)
        J00 = np.einsum('i,ni->n', dxi, xy[:, :, 0])
        J01 = np.einsum('i,ni->n', dxi, xy[:, :, 1])
        J10 = np.einsum('i,ni->n', deta, xy[:, :, 0])
        J11 = np.einsum('i,ni->n', deta, xy[:, :, 1])
        det = J00*J11 - J01*J10
        inv = 1.0/det
        Nx = inv[:, None]*( J11[:, None]*dxi[None, :] - J01[:, None]*deta[None, :])
        Ny = inv[:, None]*(-J10[:, None]*dxi[None, :] + J00[:, None]*deta[None, :])
        B = np.zeros((n, 3, 16))
        B[:, 0, 0::2] = Nx
        B[:, 1, 1::2] = Ny
        B[:, 2, 0::2] = Ny
        B[:, 2, 1::2] = Nx
        DB = np.einsum('nij,njk->nik', D_e, B)
        Ke += (w*det)[:, None, None]*np.einsum('ndi,ndj->nij', B, DB)
    dof = np.empty((n, 8, 2), dtype=np.int64)
    dof[:, :, 0] = 2*elems
    dof[:, :, 1] = 2*elems+1
    dof = dof.reshape(n, 16)
    rows = np.broadcast_to(dof[:, :, None], (n, 16, 16)).ravel()
    cols = np.broadcast_to(dof[:, None, :], (n, 16, 16)).ravel()
    K = coo_matrix((Ke.ravel(), (rows, cols)), shape=(2*len(coords),)*2)
    return K.tocsr()


def elem_stress_at(coords, elem, D, u, x, y):
    """Stress (s11,s22,s12) at physical point (x,y) inside element."""
    ec = coords[elem]
    xi = eta = 0.0
    for _ in range(50):
        N, dxi, deta = q8_shape(xi, eta)
        xe = N @ ec[:, 0]; ye = N @ ec[:, 1]
        J00 = dxi @ ec[:, 0]; J01 = dxi @ ec[:, 1]
        J10 = deta @ ec[:, 0]; J11 = deta @ ec[:, 1]
        det = J00*J11 - J01*J10
        r0 = x - xe; r1 = y - ye
        dxi_ = ( J11*r0 - J01*r1)/det
        deta_ = (-J10*r0 + J00*r1)/det
        xi += dxi_; eta += deta_
        if abs(dxi_) + abs(deta_) < 1e-14:
            break
    N, dxi, deta = q8_shape(xi, eta)
    J00 = dxi @ ec[:, 0]; J01 = dxi @ ec[:, 1]
    J10 = deta @ ec[:, 0]; J11 = deta @ ec[:, 1]
    det = J00*J11 - J01*J10
    inv = 1.0/det
    Nx = ( J11*dxi - J01*deta)*inv
    Ny = (-J10*dxi + J00*deta)*inv
    ue = np.empty(16)
    ue[0::2] = u[2*elem]
    ue[1::2] = u[2*elem+1]
    # strains
    e11 = Nx @ ue[0::2]
    e22 = Ny @ ue[1::2]
    g12 = Nx @ ue[1::2] + Ny @ ue[0::2]
    return D @ np.array([e11, e22, g12])


# ---------------------------------------------------------------- mesh
def geom_breaks(total, d0, q, cap):
    """Breakpoints 0..total with first interval d0, growth q, capped at cap."""
    if d0 >= total:
        return np.array([0.0, total])
    xs = [0.0]
    s = d0
    while True:
        nxt = xs[-1] + s
        if nxt >= total - 1e-13:
            xs.append(total)
            break
        xs.append(nxt)
        s = min(s*q, cap)
    return np.array(xs)


class BilayerCrack:
    """Structured Q8 mesh; interface nodes duplicated for x1 < 0 (crack)."""

    def __init__(self, h, H, L, d1, d2, q1=1.32, q2=1.8, cap1=1.0, cap2=0.35):
        Xp = geom_breaks(L, d1, q1, cap1)
        Xm = -Xp[::-1][:-1]
        self.X = np.concatenate([Xm, Xp])
        self.it = len(Xm)                       # index of x1 = 0
        Yu = geom_breaks(h, d2, q2, cap2)
        Yl = -geom_breaks(H, d2, q2, cap2)[::-1][:-1]
        self.Y = np.concatenate([Yl, Yu])
        self.j0 = len(Yl)                       # index of x2 = 0
        nx, ny = len(self.X), len(self.Y)

        nid = np.arange(nx*ny).reshape(nx, ny)
        hbase = nx*ny
        hnid = hbase + np.arange((nx-1)*ny).reshape(nx-1, ny)
        vbase = hbase + (nx-1)*ny
        vnid = vbase + np.arange(nx*(ny-1)).reshape(nx, ny-1)
        nn = vbase + nx*(ny-1)
        coords = np.zeros((nn, 2))
        for i in range(nx):
            for j in range(ny):
                coords[nid[i, j]] = (self.X[i], self.Y[j])
        for i in range(nx-1):
            for j in range(ny):
                coords[hnid[i, j]] = (0.5*(self.X[i]+self.X[i+1]), self.Y[j])
        for i in range(nx):
            for j in range(ny-1):
                coords[vnid[i, j]] = (self.X[i], 0.5*(self.Y[j]+self.Y[j+1]))

        # duplicate interface-row nodes with x < 0 (crack faces)
        remap = np.arange(nn)
        dup_info = []          # (base_id, dup_id)
        dup_coords = []
        nxt = nn
        for i in range(nx):
            if self.X[i] < 0:
                b = nid[i, self.j0]
                dup_coords.append((self.X[i], 0.0))
                dup_info.append((b, nxt)); remap[b] = nxt; nxt += 1
        for i in range(nx-1):
            if 0.5*(self.X[i]+self.X[i+1]) < 0:
                b = hnid[i, self.j0]
                dup_coords.append((0.5*(self.X[i]+self.X[i+1]), 0.0))
                dup_info.append((b, nxt)); remap[b] = nxt; nxt += 1
        if dup_coords:
            coords = np.vstack([coords, np.array(dup_coords)])
        self.coords = coords
        self.remap = remap
        self.nid, self.hnid, self.vnid = nid, hnid, vnid
        self.nx, self.ny = nx, ny

        # elements (upper block keeps base ids; lower block uses remap)
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
        self.elem_of_cell = np.array(elems).reshape(ny-1, nx-1, 8)

        # crack-face node pairs (upper id, lower id, x1)
        pairs = []
        for (b, d) in dup_info:
            pairs.append((b, d, coords[b, 0]))
        self.pairs = pairs

        # left-end edge list: (n1,n2,mid, y_bottom, y_top, material)
        edges = []
        for j in range(ny-1):
            n = [nid[0, j], nid[0, j+1], vnid[0, j]]
            if j < self.j0:                     # lower layer: remap top node
                n = [remap[k] for k in n]
            edges.append((n, self.Y[j], self.Y[j+1], 1 if j >= self.j0 else 2))
        self.left_edges = edges


def end_loads(mesh, sigma1):
    """Consistent nodal forces for -sigma1(x2) e_x on the left end faces.

    sigma1: dict material-> callable(x2).
    """
    f = np.zeros(2*len(mesh.coords))
    for (n, yb, yt, mat) in mesh.left_edges:
        s = sigma1[mat]
        half = 0.5*(yt-yb)
        for (t, w) in ((-np.sqrt(3/5), 5/9), (0.0, 8/9), (np.sqrt(3/5), 5/9)):
            x2 = 0.5*(yb+yt) + half*t
            s11 = s(x2)
            Na = 0.5*t*(t-1); Nb = 0.5*t*(t+1); Nm = 1-t*t
            f[2*n[0]]   += -Na*s11*half*w
            f[2*n[1]]   += -Nb*s11*half*w
            f[2*n[2]]   += -Nm*s11*half*w
    return f


def solve(mesh, D_list, f, pin=None):
    K = assemble(mesh.coords, mesh.elems, mesh.mats, D_list)
    if pin is None:
        pin = [(mesh.nid[0, 0], 'both'), (mesh.nid[0, mesh.ny-1], 'y')]
    fixed = set()
    for (node, which) in pin:
        fixed.add(2*node)
        if which in ('both', 'x'):
            pass
        if which == 'both':
            fixed.add(2*node+1)
    fixed = np.array(sorted(fixed), dtype=np.int64)
    alldof = np.arange(2*len(mesh.coords))
    free = np.setdiff1d(alldof, fixed)
    u = np.zeros(2*len(mesh.coords))
    u[free] = spsolve(K[free][:, free], f[free])
    return u


def interface_tractions(mesh, D_list, u, rs):
    """T = s22 + i s12 on the interface ahead of the tip, averaged over the
    elements just above and just below the interface."""
    Dm = {1: D_list[0], 2: D_list[1]}
    out = np.zeros(len(rs), dtype=complex)
    up = mesh.elem_of_cell[mesh.j0]      # row of cells with bottom on interface
    low = mesh.elem_of_cell[mesh.j0-1]
    X = mesh.X
    for m, r in enumerate(rs):
        i = np.searchsorted(X, r, side='right') - 1
        i = min(max(i, mesh.it), len(X)-2)
        eu = up[i]; el = low[i]
        su = elem_stress_at(mesh.coords, eu, Dm[1], u, r, 0.0)
        sl = elem_stress_at(mesh.coords, el, Dm[2], u, r, 0.0)
        s22 = 0.5*(su[1]+sl[1]); s12 = 0.5*(su[2]+sl[2])
        out[m] = s22 + 1j*s12
    return out


def crack_opening(mesh, u):
    """(r, delta) with delta = du_upper - du_lower (complex) behind tip."""
    rs, ds = [], []
    for (up_id, lo_id, x1) in mesh.pairs:
        dx = u[2*up_id]-u[2*lo_id]
        dy = u[2*up_id+1]-u[2*lo_id+1]
        rs.append(-x1); ds.append(dx + 1j*dy)
    o = np.argsort(rs)
    return np.array(rs)[o], np.array(ds)[o]
