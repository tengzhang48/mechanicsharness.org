# Visible messages and commands of codex_nodef1, in order

Cut from the log of the Codex CLI (codex.log, kept outside git). Commands are shortened; their output is left out.


## Message 1

I’ll read `PUBLIC_TASK.md`, work through the problem in this folder, and check what I can without using the web.

- `/usr/bin/bash -lc "pwd && rg --files -g 'PUBLIC_TASK.md' -g 'AGENTS.md' -g '"'!node_modules'"' -g '"'!vendor'"' -g '"'!'"\\.git'" in /work` (succeeded in 0ms:)
- `/usr/bin/bash -lc 'cat PUBLIC_TASK.md' in /work` (succeeded in 0ms:)
- `/usr/bin/bash -lc 'ls -la' in /work` (succeeded in 0ms:)
- `/usr/bin/bash -lc "python - <<'PY' import sys, importlib.util import numpy, scipy print('Python:',sys.version) print('NumPy:',numpy.__version__,'SciPy:',scipy.__version__) for name in ['skfem','meshio','gmsh','sympy','matplotlib']: print(name, bool(importlib.util.find_spec(name))) PY" in /work` (exited 127 in 0ms:)
- `/usr/bin/bash -lc "python3 - <<'PY' import sys, importlib.util import numpy, scipy print('Python:',sys.version) print('NumPy:',numpy.__version__,'SciPy:',scipy.__version__) for name in ['skfem','meshio','gmsh','sympy','matplotlib']: print(name, bool(importlib.util.find_spec(name))) PY" in /work` (succeeded in 40ms:)

## Message 2

The far-field beam energies give the energy release rates directly. The phase angles need a local elasticity solution, so I’m building a finite-element calculation and will check its convergence near the crack tip.
**Preparing Q2 derivative matrices**
**Deriving Q2 stiffness blocks**
**Designing a graded mesh**
**Assessing radial mesh options**
**Defining radial scaling**
**Planning radial mesh**
**Defining radial mesh topology**
**Estimating mesh resolution**
**Integrating left-edge loads**
**Selecting right-edge constraints**
**Building shared angular nodes**
**Building collapsed core mesh**
**Defining polar radial layers**
**Defining adaptive angular breaks**
**Refining top-edge angular spacing**
**Assessing Q2 far-field accuracy**
**Building the sparse solver**
**Fitting K over radii**
**Checking scaled traction fit**
**Creating solve.py scaffold**
apply patch
patch: completed
/work/solve.py
diff --git a/solve.py b/solve.py
new file mode 100644
index 0000000000000000000000000000000000000000..05d07727529a1a8bbf23391ea66b01016fce7ed1
--- /dev/null
+++ b/solve.py
@@ -0,0 +1,254 @@
+"""Plane-strain Q2 finite elements for the traction-free interface crack.
+
+Only NumPy/SciPy are used.  All lengths are in h, and E2=P=M=1 in the
+separate load cases.  The mesh consists of graded radial rings with its
+outer ring on the rectangle [-L,L] x [-2,1].
+"""
+
+import argparse
+import json
+import time
+from pathlib import Path
+
+import numpy as np
+from numpy.polynomial.legendre import leggauss
+from scipy.sparse import coo_matrix
+from scipy.sparse.linalg import splu
+
+
+NU = 1.0 / 3.0
+H1, H2 = 1.0, 2.0
+BETA = 0.2
+EPSILON = np.log((1 - BETA) / (1 + BETA)) / (2 * np.pi)
+G_COEFFICIENT = (1 - BETA**2) * (1 - NU**2) * (1 + 1 / 9) / 2
+
+
+def basis(s):
+    s = np.asarray(s)
+    n = np.stack((s * (s - 1) / 2, 1 - s * s, s * (s + 1) / 2), axis=-1)
+    d = np.stack((s - 0.5, -2 * s, s + 0.5), axis=-1)
+    return n, d
+
+
+def shape(xi, eta):
+    a, da = basis(xi)
+    b, db = basis(eta)
+    n = (b[..., :, None] * a[..., None, :]).reshape(-1, 9)
+    dx = (b[..., :, None] * da[..., None, :]).reshape(-1, 9)
+    dy = (db[..., :, None] * a[..., None, :]).reshape(-1, 9)
+    return n, dx, dy
+
+
+def mesh(length, angle_deg, ratio, rmin):
+    # Split at corners and at the two material interfaces, so every outer
+    # element edge is straight and every element belongs to one material.
+    a = np.arctan(H1 / length)
+    b = np.arctan(H2 / length)
+    corners = [-np.pi, -np.pi + b, -np.pi / 2, -b, 0,
+               a, np.pi / 2, np.pi - a, np.pi]
+    coarse_theta = []
+    for start, stop in zip(corners[:-1], corners[1:]):
+        count = max(1, int(np.ceil((stop - start) / np.deg2rad(angle_deg))))
+        coarse_theta.extend(np.linspace(start, stop, count + 1)[:-1])
+    coarse_theta.append(np.pi)
+    coarse_theta = np.array(coarse_theta)
+    theta = np.empty(2 * len(coarse_theta) - 1)
+    theta[::2] = coarse_theta
+    theta[1::2] = (coarse_theta[:-1] + coarse_theta[1:]) / 2
+    # All interior rings are circular; beyond radius .4 they smoothly turn
+    # into the rectangular outer boundary.
+    r_inner = 0.4
+    count = int(np.ceil(np.log(r_inner / rmin) / np.log(ratio)))
+    rcoarse = np.r_[0, np.geomspace(rmin, r_inner, count + 1)]
+    rr = np.empty(2 * len(rcoarse) - 1)
+    rr[::2] = rcoarse
+    rr[1::2] = (rcoarse[:-1] + rcoarse[1:]) / 2
+    out_count = int(np.ceil(np.log(np.hypot(length, H2) / r_inner) / np.log(ratio)))
+    tout = np.linspace(0, 1, 2 * out_count + 1)[1:]
+    ct, st = np.cos(theta), np.sin(theta)
+    outer = np.minimum(length / np.maximum(np.abs(ct), 1e-100),
+                       np.where(st >= 0, H1, H2) / np.maximum(np.abs(st), 1e-100))
+    radii = np.concatenate((np.broadcast_to(rr[:, None], (len(rr), len(theta))),
+                            r_inner * (outer[None, :] / r_inner)**tout[:, None]), axis=0)
+    xx, yy = radii * ct, radii * st
+    # Remove floating point noise on material boundaries and rectangle sides.
+    yy[np.abs(yy) < 1e-14 * np.maximum(radii, 1e-100)] = 0
+    xx[np.abs(xx) < 1e-14 * np.maximum(radii, 1e-100)] = 0
+    nodes = np.column_stack((xx[1:].ravel(), yy[1:].ravel()))
+    nodes = np.vstack(([0, 0], nodes))
+    ids = np.zeros(xx.shape, dtype=np.int64)
+    ids[1:] = np.arange(1, len(nodes)).reshape(xx.shape[0] - 1, xx.shape[1])
+    elements, ring_ids, angle_ids = [], [], []
+    for i in range((len(radii) - 1) // 2):
+        for j in range((len(theta) - 1) // 2):
+            # Local ordering: angle first, radius second.
+            elements.append(ids[2 * i:2 * i + 3, 2 * j:2 * j + 3].T.ravel())
+            ring_ids.append(i)
+            angle_ids.append(j)
+    elements = np.array(elements)
+    ring_ids, angle_ids = np.array(ring_ids), np.array(angle_ids)
+    modulus = np.where(coarse_theta[angle_ids] >= -1e-14, 9., 1.)
+    return dict(nodes=nodes, elements=elements, modulus=modulus, ids=ids,
+                radii=radii, theta=theta, coarse_theta=coarse_theta,
+                ring_ids=ring_ids, angle_ids=angle_ids, nr_inner=len(rcoarse) - 1)
+
+
+def gradients(coords, dxi, deta):
+    xxi = np.einsum('qa,eac->eqc', dxi, coords)
+    xeta = np.einsum('qa,eac->eqc', deta, coords)
+    det = xxi[:, :, 0] * xeta[:, :, 1] - xxi[:, :, 1] * xeta[:, :, 0]
+    dx = (xeta[:, :, 1, None] * dxi[None, :, :] -
+          xxi[:, :, 1, None] * deta[None, :, :]) / det[:, :, None]
+    dy = (-xeta[:, :, 0, None] * dxi[None, :, :] +
+          xxi[:, :, 0, None] * deta[None, :, :]) / det[:, :, None]
+    return dx, dy, det
+
+
+def assemble(m, quadrature=4):
+    nodes, elements = m['nodes'], m['elements']
+    gp, gw = leggauss(quadrature)
+    xi, eta = np.meshgrid(gp, gp)
+    n, dxi, deta = shape(xi.ravel(), eta.ravel())
+    weights = (gw[:, None] * gw[None, :]).ravel()
+    ne = len(elements)
+    ke = np.zeros((ne, 18, 18))
+    min_det = np.inf
+    for begin in range(0, ne, 1500):
+        stop = min(begin + 1500, ne)
+        coords = nodes[elements[begin:stop]]
+        dx, dy, det = gradients(coords, dxi, deta)
+        min_det = min(min_det, det.min())
+        if np.any(det <= 0):
+            raise ValueError('Inverted element')
+        emod = m['modulus'][begin:stop]
+        mu = emod / (2 * (1 + NU))
+        lam = emod * NU / ((1 + NU) * (1 - 2 * NU))
+        wt = weights[None, :] * det
+        ddx = np.einsum('eqa,eqb,eq->eab', dx, dx, wt, optimize=True)
+        ddy = np.einsum('eqa,eqb,eq->eab', dy, dy, wt, optimize=True)
+        dxy = np.einsum('eqa,eqb,eq->eab', dx, dy, wt, optimize=True)
+        k = ke[begin:stop]
+        k[:, ::2, ::2] = (lam + 2 * mu)[:, None, None] * ddx + mu[:, None, None] * ddy
+        k[:, 1::2, 1::2] = mu[:, None, None] * ddx + (lam + 2 * mu)[:, None, None] * ddy
+        k[:, ::2, 1::2] = lam[:, None, None] * dxy + mu[:, None, None] * dxy.transpose(0, 2, 1)
+        k[:, 1::2, ::2] = k[:, ::2, 1::2].transpose(0, 2, 1)
+    dofs = np.stack((2 * elements, 2 * elements + 1), axis=-1).reshape(ne, 18)
+    ii = np.broadcast_to(dofs[:, :, None], ke.shape).ravel()
+    jj = np.broadcast_to(dofs[:, None, :], ke.shape).ravel()
+    kglobal = coo_matrix((ke.ravel(), (ii, jj)), shape=(2 * len(nodes), 2 * len(nodes))).tocsc()
+
+    f = np.zeros((2 * len(nodes), 2))
+    # Prescribe the stated axial traction on the remote left edge.
+    for j in range((len(m['theta']) - 1) // 2):
+        edge = m['ids'][-1, 2 * j:2 * j + 3]
+        xyz = nodes[edge]
+        if np.all(np.abs(xyz[:, 0] + np.max(nodes[:, 0])) < 1e-8):
+            nn, dn = basis(gp)
+            y = nn @ xyz[:, 1]
+            jac = np.linalg.norm(dn @ xyz, axis=1)
+            if y.mean() > 0:
+                traction = np.column_stack((np.ones_like(y), 12 * (y - .5)))
+            else:
+                traction = np.column_stack((-.5 - 12 * 1.5 * (y + 1) / 8,
+                                             -12 * (y + 1) / 8))
+            fe = nn.T @ ((gw * jac)[:, None] * traction)
+            f[2 * edge] += fe
+    # The bonded material far ahead is stress-free. Fix its remote edge.
+    right = np.flatnonzero(np.abs(nodes[:, 0] - np.max(nodes[:, 0])) < 1e-8)
+    constrained = np.ravel(np.column_stack((2 * right, 2 * right + 1)))
+    free = np.setdiff1d(np.arange(2 * len(nodes)), constrained)
+    return kglobal, f, free, min_det
+
+
+def tractions(m, u, side):
+    # Use three sample points per radial element on the bonded interface.
+    # Discard the innermost core and all outer, non-circular rings.
+    j0 = np.argmin(np.abs(m['coarse_theta']))
+    j = j0 if side == 'upper' else j0 - 1
+    selected = np.flatnonzero((m['angle_ids'] == j) & (m['ring_ids'] > 0) &
+                              (m['ring_ids'] < m['nr_inner']))
+    xi = np.array([-.6, 0., .6])
+    eta = np.full_like(xi, -1. if side == 'upper' else 1.)
+    n, dxi, deta = shape(xi, eta)
+    el = m['elements'][selected]
+    coords = m['nodes'][el]
+    dx, dy, det = gradients(coords, dxi, deta)
+    ue = u[np.stack((2 * el, 2 * el + 1), axis=-1)]
+    ux = np.einsum('eqa,eack->eqck', dx, ue)
+    uy = np.einsum('eqa,eack->eqck', dy, ue)
+    e = m['modulus'][selected, None, None]
+    mu, lam = e / (2 * (1 + NU)), e * NU / ((1 + NU) * (1 - 2 * NU))
+    sig22 = lam * ux[:, :, 0] + (lam + 2 * mu) * uy[:, :, 1]
+    sig12 = mu * (uy[:, :, 0] + ux[:, :, 1])
+    t = sig22 + 1j * sig12
+    r = np.einsum('qa,eac->eqc', n, coords)[:, :, 0]
+    return r.ravel(), t.reshape(-1, 2)
+
+
+def exact_values():
+    ga = (1 - NU**2) * (.5 / 9 + .5 / 2 + 6 * 1.5**2 / 2**3)
+    gb = (1 - NU**2) * (6 / 9 + 6 / 2**3)
+    gab = (1 - NU**2) * 6 * 1.5 / 2**3
+    return ga, gb, gab
+
+
+def fit_tip(r, t, rlo, rhi):
+    select = (r >= rlo) & (r <= rhi)
+    rr = r[select]
+    kh = t[select] * (np.sqrt(2 * np.pi * rr) * np.exp(-1j * EPSILON * np.log(rr)))[:, None]
+    design = np.column_stack((np.ones_like(rr), rr / rhi, (rr / rhi)**2))
+    fitted, *_ = np.linalg.lstsq(design, kh, rcond=None)
+    return fitted[0], np.sqrt(np.mean(np.abs(design @ fitted - kh)**2, axis=0))
+
+
+def run(args):
+    start = time.monotonic()
+    m = mesh(args.length, args.angle, args.ratio, args.rmin)
+    print(f"Mesh: {len(m['nodes'])} nodes, {len(m['elements'])} Q2 elements", flush=True)
+    k, f, free, min_det = assemble(m, args.quadrature)
+    print(f'Assembled in {time.monotonic() - start:.1f}s; solving {len(free)} unknowns', flush=True)
+    factor = splu(k[free][:, free])
+    u = np.zeros_like(f)
+    u[free] = factor.solve(f[free])
+    residual = k @ u - f
+    work = f.T @ u
+    print(f'Solved in {time.monotonic() - start:.1f}s; relative residual {np.linalg.norm(residual[free]) / np.linalg.norm(f[free]):.3e}', flush=True)
+    print('Remote force resultants:', f[::2].sum(axis=0), f[1::2].sum(axis=0), flush=True)
+    print('Remote moments:', -(m['nodes'][:, 1, None] * f[::2]).sum(axis=0), flush=True)
+    data = dict(parameters=vars(args), elapsed_seconds=time.monotonic() - start,
+                nodes=len(m['nodes']), elements=len(m['elements']), work_matrix=work.tolist(),
+                epsilon=EPSILON, exact_G=exact_values(), fits={})
+    all_r, all_t = [], []
+    for side in ['upper', 'lower']:
+        r, t = tractions(m, u, side)
+        all_r.append(r)
+        all_t.append(t)
+    # Both materials should give the same traction. Keep the separate
+    # evaluations as a discretization check as well as their average.
+    for name, r, t in [('upper', all_r[0], all_t[0]), ('lower', all_r[1], all_t[1]),
+                        ('average', all_r[0], (all_t[0] + all_t[1]) / 2)]:
+        for rlo, rhi in [(1e-4, 1e-2), (3e-5, 3e-3), (1e-5, 1e-3), (3e-6, 3e-4)]:
+            if rlo < args.rmin * 15:
+                continue
+            kh, error = fit_tip(r, t, rlo, rhi)
+            key = f'{name}_{rlo:g}_{rhi:g}'
+            data['fits'][key] = dict(amplitude=np.abs(kh).tolist(), phase=np.angle(kh, deg=True).tolist(),
+                                     both_phase=float(np.angle(kh.sum(), deg=True)),
+                                     residual=error.tolist(), real=kh.real.tolist(), imag=kh.imag.tolist())
+            print(key, 'K=', np.abs(kh), 'phase=', np.angle(kh, deg=True),
+                  'both=', np.angle(kh.sum(), deg=True), flush=True)
+    outfile = Path(args.output)
+    outfile.write_text(json.dumps(data, indent=2) + '\n')
+    np.savez(outfile.with_suffix('.npz'), r=all_r[0], upper=all_t[0], lower=all_t[1])
+    print(f'Wrote {outfile}; total {time.monotonic() - start:.1f}s', flush=True)
+
+
+if __name__ == '__main__':
+    parser = argparse.ArgumentParser(description=__doc__)
+    parser.add_argument('--length', type=float, default=12)
+    parser.add_argument('--angle', type=float, default=5)
+    parser.add_argument('--ratio', type=float, default=1.22)
+    parser.add_argument('--rmin', type=float, default=1e-7)
+    parser.add_argument('--quadrature', type=int, default=4)
+    parser.add_argument('--output', default='fem_initial.json')
+    run(parser.parse_args())

**Running mesh ordering checks**
**Checking traction quadrature**

- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 solve.py --output fem_initial.json' in /work` (succeeded in 4648ms:)
- `/usr/bin/bash -lc "python3 - <<'PY' import numpy as np import solve x,y,c=solve.exact_values() print('G A,B,cross',x,y,c) print('Amplitudes',np.sqrt(np.array([x,y])/solve.G_COEFFICIENT)) print('relative phase',np.degrees(np.arccos(c/np.sqrt(x*y)))) print('phase per decade',solve.EPSILON*np.log(10)*1 ...` (succeeded in 34ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 solve.py --angle 3 --ratio 1.15 --rmin 1e-8 --output fem_refined.json' in /work **Refining core convergence**` (succeeded in 18040ms:)
- `/usr/bin/bash -lc "python3 - <<'PY' import numpy as np from scipy.optimize import least_squares from solve import EPSILON z=np.load('fem_initial.npz') r=z['r']; t=(z['upper']+z['lower'])/2 choose=(r>1e-4)&(r<1e-2) for i in range(2): x=np.log(r[choose]); y=np.unwrap(np.angle(t[choose,i])) design=np.c ...` (succeeded in 192ms:)

## Message 3

The first two meshes give phases near 52.25° for case A and 4.3° for case B. They also recover the predicted change of −8.5136° per decade in distance. I’m now checking mesh refinement and the effect of truncating the infinite layers.

- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 solve.py --length 8 --angle 3 --ratio 1.15 --rmin 1e-8 --output fem_length8.json' in /work` (succeeded in 17050ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 solve.py --length 10 --angle 2 --ratio 1.12 --rmin 1e-9 --output fem_fine.json' in /work **Preparing verification checks** **Checking crack-tip stress terms**` (succeeded in 49091ms:)

## Message 4

The finer mesh gives 52.249° for case A and 4.287° for case B. The mesh and domain-size changes shift the refined results by less than 0.03°, and the computed amplitudes agree with the energy calculation within 0.03%. I’m writing the requested files and documenting the checks.
apply patch
patch: completed
/work/postprocess.py
diff --git a/postprocess.py b/postprocess.py
new file mode 100644
index 0000000000000000000000000000000000000000..00c301ff6de463026e6bd0fbe88e62620334fb43
--- /dev/null
+++ b/postprocess.py
@@ -0,0 +1,138 @@
+"""Check the saved elasticity runs and write the requested results.json."""
+
+import json
+from pathlib import Path
+
+import numpy as np
+from numpy.polynomial.legendre import leggauss
+
+from solve import BETA, EPSILON, G_COEFFICIENT, NU, exact_values
+
+
+def read_run(name):
+    return json.loads(Path(name).read_text())
+
+
+def check_beam_energy(p, moment):
+    """Independently integrate the stated remote stress energy densities."""
+    q, w = leggauss(8)
+    y1 = (q + 1) / 2
+    y2 = q - 1
+    s1 = -p - 12 * moment * (y1 - .5)
+    s2 = p / 2 + 12 * (moment + 1.5 * p) * (y2 + 1) / 8
+    return float((1 - NU**2) / 2 * (np.dot(w / 2, s1**2) / 9 + np.dot(w, s2**2)))
+
+
+def main():
+    ga, gb, gab = exact_values()
+    amplitudes = np.sqrt(np.array([ga, gb]) / G_COEFFICIENT)
+    fine = read_run('fem_fine.json')
+    reference_key = 'average_0.0001_0.01'
+    fitted = fine['fits'][reference_key]
+    phase = np.asarray(fitted['phase'])
+    # Use the exact energy-based amplitudes in superposing the two unit loads.
+    # h=1, so M=P h corresponds to adding these two complex coefficients.
+    kh = amplitudes * np.exp(1j * np.deg2rad(phase))
+    both_phase = float(np.angle(kh.sum(), deg=True))
+
+    results = {
+        'caseA': {
+            'G_E2_h_over_P2': ga,
+            'amplitude_sqrt_h_over_P': float(amplitudes[0]),
+            'phase_at_h_deg': round(float(phase[0]), 2),
+        },
+        'caseB': {
+            'G_E2_h3_over_M2': gb,
+            'amplitude_h32_over_M': float(amplitudes[1]),
+            'phase_at_h_deg': round(float(phase[1]), 2),
+        },
+        'phase_at_h_both_deg': round(both_phase, 2),
+        'phase_change_per_decade_deg': float(np.rad2deg(EPSILON * np.log(10))),
+    }
+
+    integrated_energies = [check_beam_energy(1, 0), check_beam_energy(0, 1),
+                           check_beam_energy(1, 1)]
+    analytic_energies = [ga, gb, ga + gb + 2 * gab]
+    assert np.allclose(integrated_energies, analytic_energies, rtol=1e-13)
+
+    # Polarization of G(P,M) determines the cosine of the phase difference.
+    relative_phase_exact = float(np.rad2deg(np.arccos(gab / np.sqrt(ga * gb))))
+    relative_phase_fem = float(phase[0] - phase[1])
+    relative_phase_error = relative_phase_fem - relative_phase_exact
+    assert abs(relative_phase_error) < .05
+
+    amplitude_error = (np.asarray(fitted['amplitude']) / amplitudes - 1) * 100
+    assert np.max(np.abs(amplitude_error)) < .1
+
+    upper = np.asarray(fine['fits']['upper_0.0001_0.01']['phase'])
+    lower = np.asarray(fine['fits']['lower_0.0001_0.01']['phase'])
+    assert np.max(np.abs(upper - lower)) < .005
+
+    refined = read_run('fem_refined.json')
+    shorter = read_run('fem_length8.json')
+    initial = read_run('fem_initial.json')
+    refinement_change = phase - np.asarray(refined['fits'][reference_key]['phase'])
+    domain_change = (np.asarray(refined['fits'][reference_key]['phase']) -
+                     np.asarray(shorter['fits'][reference_key]['phase']))
+    assert np.max(np.abs(refinement_change)) < .03
+    assert np.max(np.abs(domain_change)) < .02
+
+    fit_window_change = (np.asarray(fine['fits']['average_1e-05_0.001']['phase']) - phase)
+    assert np.max(np.abs(fit_window_change)) < .005
+
+    # Estimate the oscillation from the computed stresses without imposing
+    # epsilon. Terms proportional to r and r^2 account for regular corrections.
+    z = np.load('fem_fine.npz')
+    r, t = z['r'], (z['upper'] + z['lower']) / 2
+    selected = (r > 1e-4) & (r < 1e-2)
+    rr = r[selected]
+    design = np.column_stack((np.ones_like(rr), np.log(rr), rr, rr**2))
+    fitted_epsilon = []
+    for case in range(2):
+        unwrapped_phase = np.unwrap(np.angle(t[selected, case]))
+        coefficients = np.linalg.lstsq(design, unwrapped_phase, rcond=None)[0]
+        fitted_epsilon.append(float(coefficients[1]))
+    assert np.max(np.abs(np.array(fitted_epsilon) - EPSILON)) < 2e-4
+
+    # Check linear superposition against the mixed term of the energy.
+    combined_energy = float(G_COEFFICIENT * abs(kh.sum())**2)
+    combined_energy_error = (combined_energy / analytic_energies[2] - 1) * 100
+    assert abs(combined_energy_error) < .02
+
+    checks = dict(
+        beta=BETA, epsilon=EPSILON, G_over_amplitude_squared=G_COEFFICIENT,
+        energies_from_stress_integration=integrated_energies,
+        exact_amplitudes=amplitudes.tolist(),
+        relative_phase_from_energy_deg=relative_phase_exact,
+        relative_phase_from_fem_deg=relative_phase_fem,
+        relative_phase_error_deg=relative_phase_error,
+        amplitude_error_percent=amplitude_error.tolist(),
+        upper_minus_lower_phase_deg=(upper - lower).tolist(),
+        refinement_phase_change_deg=refinement_change.tolist(),
+        length12_minus_length8_phase_deg=domain_change.tolist(),
+        fit_window_phase_change_deg=fit_window_change.tolist(),
+        epsilon_from_fem=fitted_epsilon,
+        combined_energy_error_percent=combined_energy_error,
+        unrounded_final_phases_deg=phase.tolist(),
+        unrounded_combined_phase_deg=both_phase,
+        mesh_convergence=[dict(
+            file=name, nodes=run['nodes'], elements=run['elements'],
+            length=run['parameters']['length'],
+            maximum_angle_deg=run['parameters']['angle'],
+            radial_ratio=run['parameters']['ratio'],
+            minimum_radius=run['parameters']['rmin'],
+            phase_deg=run['fits'][reference_key]['phase'],
+            both_phase_deg=run['fits'][reference_key]['both_phase'],
+            amplitudes=run['fits'][reference_key]['amplitude'],
+        ) for name, run in [('fem_initial.json', initial), ('fem_refined.json', refined),
+                            ('fem_length8.json', shorter), ('fem_fine.json', fine)]],
+        all_checks_passed=True,
+    )
+    Path('checks.json').write_text(json.dumps(checks, indent=2) + '\n')
+    Path('results.json').write_text(json.dumps(results, indent=2) + '\n')
+    print(json.dumps(results, indent=2))
+    print('All checks passed; details in checks.json.')
+
+
+if __name__ == '__main__':
+    main()
diff --git a/solve.py b/solve.py
new file mode 100644
index 0000000000000000000000000000000000000000..05d07727529a1a8bbf23391ea66b01016fce7ed1
--- /dev/null
+++ b/solve.py
@@ -0,0 +1,254 @@
+"""Plane-strain Q2 finite elements for the traction-free interface crack.
+
+Only NumPy/SciPy are used.  All lengths are in h, and E2=P=M=1 in the
+separate load cases.  The mesh consists of graded radial rings with its
+outer ring on the rectangle [-L,L] x [-2,1].
+"""
+
+import argparse
+import json
+import time
+from pathlib import Path
+
+import numpy as np
+from numpy.polynomial.legendre import leggauss
+from scipy.sparse import coo_matrix
+from scipy.sparse.linalg import splu
+
+
+NU = 1.0 / 3.0
+H1, H2 = 1.0, 2.0
+BETA = 0.2
+EPSILON = np.log((1 - BETA) / (1 + BETA)) / (2 * np.pi)
+G_COEFFICIENT = (1 - BETA**2) * (1 - NU**2) * (1 + 1 / 9) / 2
+
+
+def basis(s):
+    s = np.asarray(s)
+    n = np.stack((s * (s - 1) / 2, 1 - s * s, s * (s + 1) / 2), axis=-1)
+    d = np.stack((s - 0.5, -2 * s, s + 0.5), axis=-1)
+    return n, d
+
+
+def shape(xi, eta):
+    a, da = basis(xi)
+    b, db = basis(eta)
+    n = (b[..., :, None] * a[..., None, :]).reshape(-1, 9)
+    dx = (b[..., :, None] * da[..., None, :]).reshape(-1, 9)
+    dy = (db[..., :, None] * a[..., None, :]).reshape(-1, 9)
+    return n, dx, dy
+
+
+def mesh(length, angle_deg, ratio, rmin):
+    # Split at corners and at the two material interfaces, so every outer
+    # element edge is straight and every element belongs to one material.
+    a = np.arctan(H1 / length)
+    b = np.arctan(H2 / length)
+    corners = [-np.pi, -np.pi + b, -np.pi / 2, -b, 0,
+               a, np.pi / 2, np.pi - a, np.pi]
+    coarse_theta = []
+    for start, stop in zip(corners[:-1], corners[1:]):
+        count = max(1, int(np.ceil((stop - start) / np.deg2rad(angle_deg))))
+        coarse_theta.extend(np.linspace(start, stop, count + 1)[:-1])
+    coarse_theta.append(np.pi)
+    coarse_theta = np.array(coarse_theta)
+    theta = np.empty(2 * len(coarse_theta) - 1)
+    theta[::2] = coarse_theta
+    theta[1::2] = (coarse_theta[:-1] + coarse_theta[1:]) / 2
+    # All interior rings are circular; beyond radius .4 they smoothly turn
+    # into the rectangular outer boundary.
+    r_inner = 0.4
+    count = int(np.ceil(np.log(r_inner / rmin) / np.log(ratio)))
+    rcoarse = np.r_[0, np.geomspace(rmin, r_inner, count + 1)]
+    rr = np.empty(2 * len(rcoarse) - 1)
+    rr[::2] = rcoarse
+    rr[1::2] = (rcoarse[:-1] + rcoarse[1:]) / 2
+    out_count = int(np.ceil(np.log(np.hypot(length, H2) / r_inner) / np.log(ratio)))
+    tout = np.linspace(0, 1, 2 * out_count + 1)[1:]
+    ct, st = np.cos(theta), np.sin(theta)
+    outer = np.minimum(length / np.maximum(np.abs(ct), 1e-100),
+                       np.where(st >= 0, H1, H2) / np.maximum(np.abs(st), 1e-100))
+    radii = np.concatenate((np.broadcast_to(rr[:, None], (len(rr), len(theta))),
+                            r_inner * (outer[None, :] / r_inner)**tout[:, None]), axis=0)
+    xx, yy = radii * ct, radii * st
+    # Remove floating point noise on material boundaries and rectangle sides.
+    yy[np.abs(yy) < 1e-14 * np.maximum(radii, 1e-100)] = 0
+    xx[np.abs(xx) < 1e-14 * np.maximum(radii, 1e-100)] = 0
+    nodes = np.column_stack((xx[1:].ravel(), yy[1:].ravel()))
+    nodes = np.vstack(([0, 0], nodes))
+    ids = np.zeros(xx.shape, dtype=np.int64)
+    ids[1:] = np.arange(1, len(nodes)).reshape(xx.shape[0] - 1, xx.shape[1])
+    elements, ring_ids, angle_ids = [], [], []
+    for i in range((len(radii) - 1) // 2):
+        for j in range((len(theta) - 1) // 2):
+            # Local ordering: angle first, radius second.
+            elements.append(ids[2 * i:2 * i + 3, 2 * j:2 * j + 3].T.ravel())
+            ring_ids.append(i)
+            angle_ids.append(j)
+    elements = np.array(elements)
+    ring_ids, angle_ids = np.array(ring_ids), np.array(angle_ids)
+    modulus = np.where(coarse_theta[angle_ids] >= -1e-14, 9., 1.)
+    return dict(nodes=nodes, elements=elements, modulus=modulus, ids=ids,
+                radii=radii, theta=theta, coarse_theta=coarse_theta,
+                ring_ids=ring_ids, angle_ids=angle_ids, nr_inner=len(rcoarse) - 1)
+
+
+def gradients(coords, dxi, deta):
+    xxi = np.einsum('qa,eac->eqc', dxi, coords)
+    xeta = np.einsum('qa,eac->eqc', deta, coords)
+    det = xxi[:, :, 0] * xeta[:, :, 1] - xxi[:, :, 1] * xeta[:, :, 0]
+    dx = (xeta[:, :, 1, None] * dxi[None, :, :] -
+          xxi[:, :, 1, None] * deta[None, :, :]) / det[:, :, None]
+    dy = (-xeta[:, :, 0, None] * dxi[None, :, :] +
+          xxi[:, :, 0, None] * deta[None, :, :]) / det[:, :, None]
+    return dx, dy, det
+
+
+def assemble(m, quadrature=4):
+    nodes, elements = m['nodes'], m['elements']
+    gp, gw = leggauss(quadrature)
+    xi, eta = np.meshgrid(gp, gp)
+    n, dxi, deta = shape(xi.ravel(), eta.ravel())
+    weights = (gw[:, None] * gw[None, :]).ravel()
+    ne = len(elements)
+    ke = np.zeros((ne, 18, 18))
+    min_det = np.inf
+    for begin in range(0, ne, 1500):
+        stop = min(begin + 1500, ne)
+        coords = nodes[elements[begin:stop]]
+        dx, dy, det = gradients(coords, dxi, deta)
+        min_det = min(min_det, det.min())
+        if np.any(det <= 0):
+            raise ValueError('Inverted element')
+        emod = m['modulus'][begin:stop]
+        mu = emod / (2 * (1 + NU))
+        lam = emod * NU / ((1 + NU) * (1 - 2 * NU))
+        wt = weights[None, :] * det
+        ddx = np.einsum('eqa,eqb,eq->eab', dx, dx, wt, optimize=True)
+        ddy = np.einsum('eqa,eqb,eq->eab', dy, dy, wt, optimize=True)
+        dxy = np.einsum('eqa,eqb,eq->eab', dx, dy, wt, optimize=True)
+        k = ke[begin:stop]
+        k[:, ::2, ::2] = (lam + 2 * mu)[:, None, None] * ddx + mu[:, None, None] * ddy
+        k[:, 1::2, 1::2] = mu[:, None, None] * ddx + (lam + 2 * mu)[:, None, None] * ddy
+        k[:, ::2, 1::2] = lam[:, None, None] * dxy + mu[:, None, None] * dxy.transpose(0, 2, 1)
+        k[:, 1::2, ::2] = k[:, ::2, 1::2].transpose(0, 2, 1)
+    dofs = np.stack((2 * elements, 2 * elements + 1), axis=-1).reshape(ne, 18)
+    ii = np.broadcast_to(dofs[:, :, None], ke.shape).ravel()
+    jj = np.broadcast_to(dofs[:, None, :], ke.shape).ravel()
+    kglobal = coo_matrix((ke.ravel(), (ii, jj)), shape=(2 * len(nodes), 2 * len(nodes))).tocsc()
+
+    f = np.zeros((2 * len(nodes), 2))
+    # Prescribe the stated axial traction on the remote left edge.
+    for j in range((len(m['theta']) - 1) // 2):
+        edge = m['ids'][-1, 2 * j:2 * j + 3]
+        xyz = nodes[edge]
+        if np.all(np.abs(xyz[:, 0] + np.max(nodes[:, 0])) < 1e-8):
+            nn, dn = basis(gp)
+            y = nn @ xyz[:, 1]
+            jac = np.linalg.norm(dn @ xyz, axis=1)
+            if y.mean() > 0:
+                traction = np.column_stack((np.ones_like(y), 12 * (y - .5)))
+            else:
+                traction = np.column_stack((-.5 - 12 * 1.5 * (y + 1) / 8,
+                                             -12 * (y + 1) / 8))
+            fe = nn.T @ ((gw * jac)[:, None] * traction)
+            f[2 * edge] += fe
+    # The bonded material far ahead is stress-free. Fix its remote edge.
+    right = np.flatnonzero(np.abs(nodes[:, 0] - np.max(nodes[:, 0])) < 1e-8)
+    constrained = np.ravel(np.column_stack((2 * right, 2 * right + 1)))
+    free = np.setdiff1d(np.arange(2 * len(nodes)), constrained)
+    return kglobal, f, free, min_det
+
+
+def tractions(m, u, side):
+    # Use three sample points per radial element on the bonded interface.
+    # Discard the innermost core and all outer, non-circular rings.
+    j0 = np.argmin(np.abs(m['coarse_theta']))
+    j = j0 if side == 'upper' else j0 - 1
+    selected = np.flatnonzero((m['angle_ids'] == j) & (m['ring_ids'] > 0) &
+                              (m['ring_ids'] < m['nr_inner']))
+    xi = np.array([-.6, 0., .6])
+    eta = np.full_like(xi, -1. if side == 'upper' else 1.)
+    n, dxi, deta = shape(xi, eta)
+    el = m['elements'][selected]
+    coords = m['nodes'][el]
+    dx, dy, det = gradients(coords, dxi, deta)
+    ue = u[np.stack((2 * el, 2 * el + 1), axis=-1)]
+    ux = np.einsum('eqa,eack->eqck', dx, ue)
+    uy = np.einsum('eqa,eack->eqck', dy, ue)
+    e = m['modulus'][selected, None, None]
+    mu, lam = e / (2 * (1 + NU)), e * NU / ((1 + NU) * (1 - 2 * NU))
+    sig22 = lam * ux[:, :, 0] + (lam + 2 * mu) * uy[:, :, 1]
+    sig12 = mu * (uy[:, :, 0] + ux[:, :, 1])
+    t = sig22 + 1j * sig12
+    r = np.einsum('qa,eac->eqc', n, coords)[:, :, 0]
+    return r.ravel(), t.reshape(-1, 2)
+
+
+def exact_values():
+    ga = (1 - NU**2) * (.5 / 9 + .5 / 2 + 6 * 1.5**2 / 2**3)
+    gb = (1 - NU**2) * (6 / 9 + 6 / 2**3)
+    gab = (1 - NU**2) * 6 * 1.5 / 2**3
+    return ga, gb, gab
+
+
+def fit_tip(r, t, rlo, rhi):
+    select = (r >= rlo) & (r <= rhi)
+    rr = r[select]
+    kh = t[select] * (np.sqrt(2 * np.pi * rr) * np.exp(-1j * EPSILON * np.log(rr)))[:, None]
+    design = np.column_stack((np.ones_like(rr), rr / rhi, (rr / rhi)**2))
+    fitted, *_ = np.linalg.lstsq(design, kh, rcond=None)
+    return fitted[0], np.sqrt(np.mean(np.abs(design @ fitted - kh)**2, axis=0))
+
+
+def run(args):
+    start = time.monotonic()
+    m = mesh(args.length, args.angle, args.ratio, args.rmin)
+    print(f"Mesh: {len(m['nodes'])} nodes, {len(m['elements'])} Q2 elements", flush=True)
+    k, f, free, min_det = assemble(m, args.quadrature)
+    print(f'Assembled in {time.monotonic() - start:.1f}s; solving {len(free)} unknowns', flush=True)
+    factor = splu(k[free][:, free])
+    u = np.zeros_like(f)
+    u[free] = factor.solve(f[free])
+    residual = k @ u - f
+    work = f.T @ u
+    print(f'Solved in {time.monotonic() - start:.1f}s; relative residual {np.linalg.norm(residual[free]) / np.linalg.norm(f[free]):.3e}', flush=True)
+    print('Remote force resultants:', f[::2].sum(axis=0), f[1::2].sum(axis=0), flush=True)
+    print('Remote moments:', -(m['nodes'][:, 1, None] * f[::2]).sum(axis=0), flush=True)
+    data = dict(parameters=vars(args), elapsed_seconds=time.monotonic() - start,
+                nodes=len(m['nodes']), elements=len(m['elements']), work_matrix=work.tolist(),
+                epsilon=EPSILON, exact_G=exact_values(), fits={})
+    all_r, all_t = [], []
+    for side in ['upper', 'lower']:
+        r, t = tractions(m, u, side)
+        all_r.append(r)
+        all_t.append(t)
+    # Both materials should give the same traction. Keep the separate
+    # evaluations as a discretization check as well as their average.
+    for name, r, t in [('upper', all_r[0], all_t[0]), ('lower', all_r[1], all_t[1]),
+                        ('average', all_r[0], (all_t[0] + all_t[1]) / 2)]:
+        for rlo, rhi in [(1e-4, 1e-2), (3e-5, 3e-3), (1e-5, 1e-3), (3e-6, 3e-4)]:
+            if rlo < args.rmin * 15:
+                continue
+            kh, error = fit_tip(r, t, rlo, rhi)
+            key = f'{name}_{rlo:g}_{rhi:g}'
+            data['fits'][key] = dict(amplitude=np.abs(kh).tolist(), phase=np.angle(kh, deg=True).tolist(),
+                                     both_phase=float(np.angle(kh.sum(), deg=True)),
+                                     residual=error.tolist(), real=kh.real.tolist(), imag=kh.imag.tolist())
+            print(key, 'K=', np.abs(kh), 'phase=', np.angle(kh, deg=True),
+                  'both=', np.angle(kh.sum(), deg=True), flush=True)
+    outfile = Path(args.output)
+    outfile.write_text(json.dumps(data, indent=2) + '\n')
+    np.savez(outfile.with_suffix('.npz'), r=all_r[0], upper=all_t[0], lower=all_t[1])
+    print(f'Wrote {outfile}; total {time.monotonic() - start:.1f}s', flush=True)
+
+
+if __name__ == '__main__':
+    parser = argparse.ArgumentParser(description=__doc__)
+    parser.add_argument('--length', type=float, default=12)
+    parser.add_argument('--angle', type=float, default=5)
+    parser.add_argument('--ratio', type=float, default=1.22)
+    parser.add_argument('--rmin', type=float, default=1e-7)
+    parser.add_argument('--quadrature', type=int, default=4)
+    parser.add_argument('--output', default='fem_initial.json')
+    run(parser.parse_args())

- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 postprocess.py' in /work` (succeeded in 52ms:)
- `/usr/bin/bash -lc 'cat checks.json' in /work` (succeeded in 0ms:)
- `/usr/bin/bash -lc "python3 - <<'PY' import json from pathlib import Path r = json.loads(Path('results.json').read_text()) expected = { 'caseA': {'G_E2_h_over_P2', 'amplitude_sqrt_h_over_P', 'phase_at_h_deg'}, 'caseB': {'G_E2_h3_over_M2', 'amplitude_h32_over_M', 'phase_at_h_deg'}, } assert set(r) ==  ...` (succeeded in 0ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 postprocess.py' in /work` (succeeded in 33ms:)

## Message 5

Created [results.json](/work/results.json), [method.md](/work/method.md), and reproducible solver/check code.

Phases reached: A **52.25°**, B **4.29°**, combined **30.44°**. Energy, amplitude, and convergence checks passed.

Unverified: comparison against an independent reference and rigorous numerical error bounds.

