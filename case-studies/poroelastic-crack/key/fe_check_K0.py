"""Independent FE check of problem B, loading 1: K at the tip just after loading (log P365).
Biot poroelasticity, plane strain, quarter domain [0,W]^2 with the crack on y = 0, 0 <= x <= L = 1.
Taylor-Hood (u P2, p P1), backward Euler from rest, mesh refined ~16 times towards the tip (1, 0).
Measures, at each small time t, the intensity inside the drained tip layer (drained near-tip
opening, r << sqrt(ct)) and the undrained outer intensity (sqrt(ct) << r << L); their ratio is the
contested K(0+)/K_outer: 0.75 (Detournay & Cheng Eq. 29, heuristic) or ~0.807 (GPT, three runs).
Run: <home>/miniconda3/envs/fenicsx-env/bin/python fe_check_K0.py"""
import numpy as np, ufl, basix.ufl, json, sys, time
from mpi4py import MPI
from dolfinx import mesh, fem
from dolfinx.fem.petsc import LinearProblem

import os
G, c = 1.0, 1.0
nu, nuu, alpha = float(os.environ.get("NU", 0.2)), float(os.environ.get("NUU", 0.4)), float(os.environ.get("ALPHA", 0.89))
lam = 2 * G * nu / (1 - 2 * nu)
if nuu < 0.5:
    lamu = 2 * G * nuu / (1 - 2 * nuu)
    M = (lamu - lam) / alpha ** 2                               # Biot modulus from lambda_u = lambda + alpha^2 M
    kappa = c * alpha ** 2 * (1 - 2 * nu) ** 2 * (1 - nuu) / (2 * G * (1 - nu) * (nuu - nu))   # Detournay-Cheng Eq. 8
    c_std = kappa * M * (lam + 2 * G) / (lamu + 2 * G)          # standard formula, must equal c
else:                                                           # gel limit (P371): incompressible constituents, 1/M = 0
    M = float("inf"); kappa = c / (lam + 2 * G); c_std = kappa * (lam + 2 * G)
W, L = 20.0, 1.0
NREF = int(sys.argv[1]) if len(sys.argv) > 1 else 16
GROW = float(sys.argv[2]) if len(sys.argv) > 2 else 1.25
DT0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-7
STRIP = int(sys.argv[4]) if len(sys.argv) > 4 else 0       # levels of refinement in a strip along the crack face
TMAX = float(sys.argv[5]) if len(sys.argv) > 5 else 3e-5
TAG = f"N{NREF}_g{GROW}_dt{DT0:g}" + (os.environ.get("TAGX", "")) + (f"_S{STRIP}_T{TMAX:g}" if STRIP else "")

msh = mesh.create_rectangle(MPI.COMM_WORLD, [[0, 0], [W, W]], [40, 40], mesh.CellType.triangle)
for k in range(NREF):
    R = 2.0 * 0.5 ** k
    msh.topology.create_entities(1)
    edges = mesh.locate_entities(msh, 1, lambda x: np.hypot(x[0] - L, x[1]) < R)
    out = mesh.refine(msh, edges)
    msh = out[0] if isinstance(out, tuple) else out
for k in range(STRIP):                      # resolve the diffusion layer along the whole face (P366)
    Sk = max(0.4 * 0.62 ** k, 0.02)
    msh.topology.create_entities(1)
    edges = mesh.locate_entities(msh, 1, lambda x: (x[0] < L + 0.05) & (x[1] < Sk))
    out = mesh.refine(msh, edges)
    msh = out[0] if isinstance(out, tuple) else out
hmin = min(msh.h(2, np.arange(msh.topology.index_map(2).size_local, dtype=np.int32)))
print(f"cells {msh.topology.index_map(2).size_local}, hmin {hmin:.2e}, c check {c_std:.6f}", flush=True)

el_u = basix.ufl.element("Lagrange", msh.basix_cell(), 2, shape=(2,))
el_p = basix.ufl.element("Lagrange", msh.basix_cell(), 1)
V = fem.functionspace(msh, basix.ufl.mixed_element([el_u, el_p]))
w_old = fem.Function(V)
u, p = ufl.TrialFunctions(V); v, q = ufl.TestFunctions(V)
u0, p0 = ufl.split(w_old)
eps = lambda a: ufl.sym(ufl.grad(a))
sig_eff = lambda a: 2 * G * eps(a) + lam * ufl.tr(eps(a)) * ufl.Identity(2)
dt = fem.Constant(msh, DT0)

fdim = 1
crack = mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[1], 0) & (x[0] <= L + 1e-12))
lig = mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[1], 0) & (x[0] >= L - 1e-12))
left = mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], 0))
far = mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], W) | np.isclose(x[1], W))
tags = np.hstack([crack, lig, left, far]); vals = np.hstack([np.full(len(crack), 1), np.full(len(lig), 2), np.full(len(left), 3), np.full(len(far), 4)])
order = np.argsort(tags); mt = mesh.meshtags(msh, fdim, tags[order], vals[order].astype(np.int32))
ds = ufl.Measure("ds", domain=msh, subdomain_data=mt)
dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 4})

zeta = (lambda a, b: alpha * ufl.div(a) + b / M) if M != float("inf") else (lambda a, b: alpha * ufl.div(a))
t_crack = ufl.as_vector([0.0, 1.0])            # sigma_yy = -1 on the face (loading 1), outward normal (0,-1)
a = (ufl.inner(sig_eff(u), eps(v)) - alpha * p * ufl.div(v)) * dx + (zeta(u, p) * q + dt * kappa * ufl.inner(ufl.grad(p), ufl.grad(q))) * dx
Lf = ufl.inner(t_crack, v) * ds(1) + zeta(u0, p0) * q * dx

Vu, _ = V.sub(0).collapse(); Vuy, map_uy = V.sub(0).sub(1).collapse(); Vux, map_ux = V.sub(0).sub(0).collapse(); Vp, map_p = V.sub(1).collapse()
zero_s = fem.Function(Vuy); zero_x = fem.Function(Vux); zero_p = fem.Function(Vp)
bcs = [fem.dirichletbc(zero_s, fem.locate_dofs_topological((V.sub(0).sub(1), Vuy), fdim, lig), V.sub(0).sub(1)),
       fem.dirichletbc(zero_x, fem.locate_dofs_topological((V.sub(0).sub(0), Vux), fdim, left), V.sub(0).sub(0)),
       fem.dirichletbc(zero_p, fem.locate_dofs_topological((V.sub(1), Vp), fdim, crack), V.sub(1)),
       fem.dirichletbc(zero_p, fem.locate_dofs_topological((V.sub(1), Vp), fdim, far), V.sub(1))]
prob = LinearProblem(a, Lf, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "petsc"},
                     petsc_options_prefix="b_")

xy = Vuy.tabulate_dof_coordinates()
face = np.where(np.isclose(xy[:, 1], 0) & (xy[:, 0] < L))[0]
r_face = L - xy[face, 0]; srt = np.argsort(r_face); face, r_face = face[srt], r_face[srt]


def intensities(w, t):
    uy = w.x.array[np.asarray(map_uy)[face]]
    Kd = uy * G * np.sqrt(2 * np.pi) / (2 * (1 - nu) * np.sqrt(r_face))     # drained near-tip opening relation
    Ku = uy * G * np.sqrt(2 * np.pi) / (2 * (1 - nuu) * np.sqrt(r_face))    # undrained relation
    lt = np.sqrt(c * t)
    inner = (r_face > 5 * hmin) & (r_face < 0.05 * lt)
    outer = (r_face > 5 * lt) & (r_face < 0.03)
    return (np.median(Kd[inner]) if inner.sum() >= 3 else np.nan, inner.sum(),
            np.median(Ku[outer]) if outer.sum() >= 3 else np.nan, outer.sum(), r_face, Kd, Ku)


rows, t, profiles = [], 0.0, []
w = fem.Function(V)
for step in range(400):
    t += dt.value
    w = prob.solve()
    Kin, nin, Kout, nout, rf, Kd, Ku = intensities(w, t)
    sq = np.sqrt(np.pi * L)
    uy_all = w.x.array[np.asarray(map_uy)[face]]; xs = L - r_face; o = np.argsort(xs)
    xs_s = np.concatenate([xs[o], [L]]); uy_s = np.concatenate([uy_all[o], [0.0]])
    Vol = 4 * np.trapezoid(uy_s, xs_s)          # full crack volume: two faces x two half-cracks (quarter model)
    rows.append(dict(t=t, V=float(Vol), K_inner=float(Kin / sq), n_inner=int(nin), K_outer=float(Kout / sq), n_outer=int(nout),
                     ratio=float(Kin / Kout) if np.isfinite(Kin) and np.isfinite(Kout) else None))
    print(f"t={t:.3e} K_inner/sqrt(piL)={Kin / sq:.4f} (n={nin}) K_outer/sqrt(piL)={Kout / sq:.4f} (n={nout}) ratio={rows[-1]['ratio']}", flush=True)
    for tp in (1e-5, 1e-4, 1e-3):
        if t >= tp and not any(abs(r.get("profile_t", -1) - t) < 1e-30 for r in profiles) and all(pr["tag"] != tp for pr in profiles):
            profiles.append(dict(tag=tp, profile_t=t, r=rf.tolist(), Kd=(Kd / sq).tolist(), Ku=(Ku / sq).tolist()))
    w_old.x.array[:] = w.x.array
    dt.value = dt.value * GROW
    if t > TMAX:
        break
json.dump(dict(params=dict(G=G, nu=nu, nuu=nuu, alpha=alpha, c=c, kappa=kappa, W=W, NREF=NREF, hmin=hmin), rows=rows, profiles=profiles),
          open(f"fe_check_K0_{TAG}.json", "w"), indent=1)
np.savetxt(f"fe_check_K0_profile_{TAG}.txt", np.c_[rf, Kd, Ku], header="r  K_drained_relation  K_undrained_relation (last step)")
