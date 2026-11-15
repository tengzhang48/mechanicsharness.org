"""P2-triangle plane-strain elasticity with dolfinx 0.10, plus K extraction:
  (a) interaction (M-) integral, domain form, with the numerically derived bimaterial tip field (asym.py);
  (b) least-squares fit of the FE interface tractions ahead of the tip to
      (s22 + i s12) sqrt(2 pi r) r^(-i eps) = K + A1 r + A2 r^2   (form justified by check_eigs.py).
Post-processing is done in NumPy from grad(u_h) evaluated at reference points of each cell."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
from mpi4py import MPI
import basix
import basix.ufl
import ufl
from dolfinx import mesh, fem
import dolfinx.fem.petsc
from asym import TipField


class Solution:
    pass


def solve(pts, cells, mats, load_cases, pins, degree=2):
    """mats: {1: (E,nu) for y>0, 2: (E,nu) for y<0}.
    load_cases: list of load cases; each is a list of (facet_marker(xmid)->bool array, traction(x)->ufl 2-vector).
    pins: list of (point (x,y), component) with zero displacement.
    One assembly + one LU factorisation, one solve per load case.  Returns list of Solution."""
    cel = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
    msh = mesh.create_mesh(MPI.COMM_SELF, cells, cel, pts)
    tdim = 2
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_connectivity(0, tdim)
    gx = msh.geometry.x[:, :2]
    gdm = msh.geometry.dofmap
    cen = gx[gdm].mean(axis=1)
    cmat = np.where(cen[:, 1] > 0, 1, 2)
    Q0 = fem.functionspace(msh, ("DG", 0))
    mu = fem.Function(Q0)
    lm = fem.Function(Q0)
    for k, (E, nu) in mats.items():
        sel = np.where(cmat == k)[0]
        mu.x.array[sel] = E / (2 * (1 + nu))
        lm.x.array[sel] = E * nu / ((1 + nu) * (1 - 2 * nu))
    V = fem.functionspace(msh, ("Lagrange", degree, (2,)))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    eps = lambda w: ufl.sym(ufl.grad(w))
    sig = lambda w: 2 * mu * eps(w) + lm * ufl.tr(eps(w)) * ufl.Identity(2)
    a = ufl.inner(sig(u), eps(v)) * ufl.dx
    from petsc4py import PETSc
    bf = mesh.exterior_facet_indices(msh.topology)
    mid = mesh.compute_midpoints(msh, 1, bf)
    bcs = []
    for (p, comp) in pins:
        vert = mesh.locate_entities(msh, 0, lambda X, p=p: np.isclose(X[0], p[0], atol=1e-10) & np.isclose(X[1], p[1], atol=1e-10))
        assert len(vert) == 1, (p, vert)
        dofs = fem.locate_dofs_topological(V.sub(comp), 0, vert)
        bcs.append(fem.dirichletbc(np.float64(0.0), dofs, V.sub(comp)))
    a_form = fem.form(a)
    A = dolfinx.fem.petsc.assemble_matrix(a_form, bcs=bcs)
    A.assemble()
    ksp = PETSc.KSP().create(MPI.COMM_SELF)
    ksp.setOperators(A)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    ksp.getPC().setFactorSolverType("petsc")
    x = ufl.SpatialCoordinate(msh)
    sols = []
    for tractions in load_cases:
        tags_idx, tags_val = [], []
        for k, (marker, _) in enumerate(tractions):
            sel = bf[marker(mid)]
            tags_idx.append(sel)
            tags_val.append(np.full(len(sel), k + 1, dtype=np.int32))
        idx = np.concatenate(tags_idx)
        val = np.concatenate(tags_val)
        assert len(np.unique(idx)) == len(idx), "facet tagged twice"
        order = np.argsort(idx)
        ft = mesh.meshtags(msh, 1, idx[order].astype(np.int32), val[order])
        ds = ufl.Measure("ds", domain=msh, subdomain_data=ft, metadata={"quadrature_degree": 6})
        L = 0
        for k, (_, tfun) in enumerate(tractions):
            L = L + ufl.dot(tfun(x), v) * ds(k + 1)
        L_form = fem.form(L)
        b = dolfinx.fem.petsc.assemble_vector(L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], bcs=[bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, bcs)
        uh = fem.Function(V)
        ksp.solve(b, uh.x.petsc_vec)
        uh.x.scatter_forward()
        # residual check of the linear solve
        r = b.copy()
        A.mult(uh.x.petsc_vec, r)
        r.axpy(-1.0, b)
        S = Solution()
        S.sig, S.eps_op, S.v = sig, eps, v
        S.msh, S.uh, S.cmat, S.mats, S.V = msh, uh, cmat, mats, V
        S.gx, S.gdm = gx, gdm
        S.ndofs = V.dofmap.index_map.size_global * 2
        S.resid = r.norm() / b.norm()
        S.nload_facets = len(idx)
        sols.append(S)
    return sols


def grad_at(S, Xref, cells):
    expr = fem.Expression(ufl.grad(S.uh), Xref)
    vals = expr.eval(S.msh, cells.astype(np.int32))
    return vals.reshape(len(cells), Xref.shape[0], 2, 2)


def stress_from_grad(G, E, nu):
    mu = E / (2 * (1 + nu))
    lm = E * nu / ((1 + nu) * (1 - 2 * nu))
    e11, e22 = G[..., 0, 0], G[..., 1, 1]
    e12 = 0.5 * (G[..., 0, 1] + G[..., 1, 0])
    tr = e11 + e22
    return 2 * mu * e11 + lm * tr, 2 * mu * e22 + lm * tr, 2 * mu * e12


def _q(r, r1, r2):
    t = np.clip((r - r1) / (r2 - r1), 0, 1)
    q = 1 - (10 * t ** 3 - 15 * t ** 4 + 6 * t ** 5)
    dq = -(30 * t ** 2 - 60 * t ** 3 + 30 * t ** 4) / (r2 - r1)
    return q, dq


def interaction_integrals(S, tip, rings, qdeg=10):
    """Returns list of dicts (r1, r2, K, J) for each ring (r1, r2) around tip=(xt,0) with crack at x<xt."""
    E1, nu1 = S.mats[1]
    E2, nu2 = S.mats[2]
    tf = TipField(E1, nu1, E2, nu2)
    mu1, k1, mu2, k2 = tf.mu1, tf.k1, tf.mu2, tf.k2
    c1, c2 = (k1 + 1) / mu1, (k2 + 1) / mu2
    fac = (c1 + c2) / (8 * np.cosh(np.pi * tf.eps) ** 2)
    Xq, wq = basix.make_quadrature(basix.CellType.triangle, qdeg)
    V0 = S.gx[S.gdm[:, 0]]
    V1 = S.gx[S.gdm[:, 1]]
    V2 = S.gx[S.gdm[:, 2]]
    rv = np.stack([np.hypot(Vk[:, 0] - tip[0], Vk[:, 1] - tip[1]) for Vk in (V0, V1, V2)], axis=1)
    out = []
    rmax_all = max(r2 for _, r2 in rings)
    rmin_all = min(r1 for r1, _ in rings)
    cand = np.where((rv.min(axis=1) < rmax_all) & (rv.max(axis=1) > rmin_all))[0]
    G = grad_at(S, Xq, cand)
    v0, v1, v2 = V0[cand], V1[cand], V2[cand]
    detJ = np.abs((v1[:, 0] - v0[:, 0]) * (v2[:, 1] - v0[:, 1]) - (v2[:, 0] - v0[:, 0]) * (v1[:, 1] - v0[:, 1]))
    xq = v0[:, None, :] + (v1 - v0)[:, None, :] * Xq[None, :, 0:1] + (v2 - v0)[:, None, :] * Xq[None, :, 1:2]
    W = detJ[:, None] * wq[None, :]
    mat = S.cmat[cand]
    s11 = np.zeros(W.shape); s22 = np.zeros(W.shape); s12 = np.zeros(W.shape)
    for k in (1, 2):
        sel = mat == k
        a, b, c = stress_from_grad(G[sel], *S.mats[k])
        s11[sel], s22[sel], s12[sel] = a, b, c
    X = xq[..., 0] - tip[0]
    Y = xq[..., 1] - tip[1]
    # keep each point on the side of its cell's material (points are interior, so sign(Y) already agrees)
    Y = np.where(mat[:, None] == 1, np.maximum(Y, 1e-300), np.minimum(Y, -1e-300))
    Fc = tf.complex_field(X.ravel(), Y.ravel())
    Fc = {k: v.reshape(X.shape) for k, v in Fc.items()}
    u11, u12, u21, u22 = G[..., 0, 0], G[..., 0, 1], G[..., 1, 0], G[..., 1, 1]
    R = np.hypot(X, Y)
    for (r1, r2) in rings:
        q, dq = _q(R, r1, r2)
        q1, q2 = dq * X / R, dq * Y / R
        Wfe = 0.5 * (s11 * u11 + s22 * u22 + s12 * (u12 + u21))
        J = np.sum(W * ((s11 * u11 + s12 * u21) * q1 + (s12 * u11 + s22 * u21) * q2 - Wfe * q1))
        Ks = []
        for Ka in (1.0, 1j):
            A = {kk: np.real(Ka * vv) for kk, vv in Fc.items()}
            ea12 = 0.5 * (A["u12"] + A["u21"])
            Wint = s11 * A["u11"] + s22 * A["u22"] + 2 * s12 * ea12
            term = ((s11 * A["u11"] + s12 * A["u21"]) * q1 + (s12 * A["u11"] + s22 * A["u21"]) * q2
                    + (A["s11"] * u11 + A["s12"] * u21) * q1 + (A["s12"] * u11 + A["s22"] * u21) * q2
                    - Wint * q1)
            Ks.append(np.sum(W * term) / fac)
        out.append(dict(r1=r1, r2=r2, K=Ks[0] + 1j * Ks[1], J=J))
    S.eps = tf.eps
    S.Gfac = (c1 + c2) / (16 * np.cosh(np.pi * tf.eps) ** 2)
    return out


def interface_tractions(S, tip, rmax):
    """Nodal-averaged FE tractions (s22, s12) at vertices on the bonded interface y=0, tip[0] < x < tip[0]+rmax,
    averaged over cells on each side; returns r, t_upper, t_lower (complex s22 + i s12)."""
    gx, gdm = S.gx, S.gdm
    vy = gx[gdm][:, :, 1]
    vx = gx[gdm][:, :, 0]
    on = (np.abs(vy) < 1e-14) & (vx > tip[0] + 1e-15) & (vx < tip[0] + rmax)
    cand = np.where(on.any(axis=1))[0]
    Xv = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    G = grad_at(S, Xv, cand)
    acc = {}
    for ic, c in enumerate(cand):
        k = S.cmat[c]
        s11, s22, s12 = stress_from_grad(G[ic], *S.mats[k])
        for lv in range(3):
            if on[c, lv]:
                key = (gdm[c, lv], k)
                acc.setdefault(key, []).append(s22[lv] + 1j * s12[lv])
    nodes = sorted(set(kk[0] for kk in acc))
    r = np.array([gx[n, 0] - tip[0] for n in nodes])
    tu = np.array([np.mean(acc[(n, 1)]) for n in nodes])
    tl = np.array([np.mean(acc[(n, 2)]) for n in nodes])
    o = np.argsort(r)
    return r[o], tu[o], tl[o]


def traction_fit(r, t, eps, rlo, rhi, nterms=3):
    sel = (r >= rlo) & (r <= rhi)
    f = t[sel] * np.sqrt(2 * np.pi * r[sel]) * r[sel] ** (-1j * eps)
    A = np.stack([r[sel] ** k for k in range(nterms)], axis=1).astype(complex)
    coef, *_ = np.linalg.lstsq(A, f, rcond=None)
    return coef[0], int(sel.sum())


def consistent_tractions(S, xtip, xmax):
    """Variationally consistent interface tractions on the bonded interface y=0, x>xtip.
    Nodal forces f_I = int_{material 1} sigma(u_h):eps(phi_I) dx at interface dofs equal
    -int_Gamma (s12, s22) phi_I ds (traction exerted by layer 2 on layer 1, outward normal -e2).
    The traction field is recovered in the P2 trace space by solving the 1-D P2 mass-matrix system
    over the whole bonded interface.  Returns r, s22 + i s12 at interface dofs with r < xmax - xtip."""
    msh = S.msh
    ncell = msh.topology.index_map(2).size_local
    ct = mesh.meshtags(msh, 2, np.arange(ncell, dtype=np.int32), S.cmat.astype(np.int32))
    dxm = ufl.Measure("dx", domain=msh, subdomain_data=ct)
    f = dolfinx.fem.petsc.assemble_vector(fem.form(ufl.inner(S.sig(S.uh), S.eps_op(S.v)) * dxm(1)))
    fa = f.array.copy().reshape(-1, 2)
    X = S.V.tabulate_dof_coordinates()[:, :2]
    on = np.where((np.abs(X[:, 1]) < 1e-14) & (X[:, 0] >= xtip - 1e-15))[0]
    # keep only dofs shared by both materials (bonded): crack-face dofs at x<xtip are excluded by x>=xtip
    o = on[np.argsort(X[on, 0])]
    xs = X[o, 0]
    n = len(xs)
    assert n % 2 == 1, n
    assert np.allclose(xs[1::2], 0.5 * (xs[0:-1:2] + xs[2::2]), rtol=0, atol=1e-12 * max(1, xs[-1]))
    import scipy.sparse as sp
    import scipy.sparse.linalg as spl
    rows, cols, vals = [], [], []
    loc = np.array([[4, 2, -1], [2, 16, 2], [-1, 2, 4]]) / 30.0
    for e in range((n - 1) // 2):
        ids = [2 * e, 2 * e + 1, 2 * e + 2]
        Le = xs[2 * e + 2] - xs[2 * e]
        for a_ in range(3):
            for b_ in range(3):
                rows.append(ids[a_]); cols.append(ids[b_]); vals.append(Le * loc[a_, b_])
    Mm = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    fi = fa[o]
    t12 = -spl.spsolve(Mm, fi[:, 0])
    t22 = -spl.spsolve(Mm, fi[:, 1])
    r = xs - xtip
    sel = (r > 0) & (r < xmax - xtip)
    return r[sel], (t22 + 1j * t12)[sel]
