# Visible messages and commands of codex_1, in order

Cut from the log of the Codex CLI (codex.log, kept outside git). Commands are shortened; their output is left out.


## Message 1

I’ll read `PUBLIC_TASK.md`, then work through the problem using only files and tools in this folder, without web access.

- `/usr/bin/bash -lc "pwd && rg --files -g 'PUBLIC_TASK.md' -g 'AGENTS.md' -g '"'!node_modules'"' -g '"'!.git'"' ." in /work` (succeeded in 0ms:)
- `/usr/bin/bash -lc 'cat PUBLIC_TASK.md' in /work` (succeeded in 0ms:)
- `/usr/bin/bash -lc "ls -la && python - <<'PY' import sys import numpy, scipy print(sys.version) print('numpy',numpy.__version__,'scipy',scipy.__version__) for s in ['skfem', 'meshio','sympy','matplotlib','numba','petsc4py']: try: m=__import__(s) print(s,getattr(m,'__version__','available')) except Im ...` (exited 127 in 0ms:)
- `/usr/bin/bash -lc "python3 - <<'PY' import sys import numpy, scipy print(sys.version) print('numpy',numpy.__version__,'scipy',scipy.__version__) for s in ['skfem', 'meshio','sympy','matplotlib','numba','petsc4py']: try: m=__import__(s) print(s,getattr(m,'__version__','available')) except ImportError ...` (succeeded in 0ms:)

## Message 2

The energy release rates follow directly from the strain energy in the two separated beams. The phase angles need a crack-tip field calculation; I’m setting up a finite-element model and will check its phases under mesh refinement and against the exact energy values.

- `/usr/bin/bash -lc "find tmp -maxdepth 2 -type f -print && python3 - <<'PY' import numpy as np nu=1/3 mu1=9/(2*(1+nu)); mu2=1/(2*(1+nu)); kap=3-4*nu eps=np.log((kap/mu1+1/mu2)/(kap/mu2+1/mu1))/(2*np.pi) c=((kap+1)/mu1+(kap+1)/mu2)/(16*np.cosh(np.pi*eps)**2) GA=(1-nu**2)/2*(1/9+1/2+12*1.5**2/8) GB=(1- ...` (succeeded in 0ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 12 --ratio 1.25 --rmin 1e-7 --length 16 --output fem_coarse' in /work` (succeeded in 2553ms:)
- `/usr/bin/bash -lc "python3 - <<'PY' import numpy as np z=np.load('fem_coarse.npz'); c=z['cod']; t=z['tractions']; j=z['j'] for r in [1e-5,1e-4,1e-3,.01,.1]: k=np.argmin(np.abs(t[:,0].real-r)); print('r=',t[k,0].real,'xi=',t[k,1].real) for label,v in [('upper',t[k,2:4]),('lower',t[k,4:6])]: print(lab ...` (succeeded in 0ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 20 --ratio 1.15 --rmin 1e-8 --length 24 --output fem_medium' in /work` (succeeded in 12283ms:)

## Message 3

The first two meshes give phases near 52.25° for case A and 4.28° for case B. Their phase difference agrees with the exact energy relation, and the refined crack-opening calculation matches the predicted stress-intensity magnitudes within 0.001%. I’m checking a finer mesh and the effect of the finite domain before writing the results.

- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 28 --ratio 1.10 --rmin 1e-8 --length 24 --quadrature 5 --output fem_fine' in /work **Improving analytic phase generation**` (succeeded in 31470ms:)
- `/usr/bin/bash -lc "ls -lh fem_*.npz fem_*.json && python3 - <<'PY' import numpy as np from solve import exact_values z=np.load('fem_medium.npz') c=z['cod'];t=z['tractions'];j=z['j'];exact=exact_values() s=(t[:,0].real>1e-4)&(t[:,0].real<1e-2) for xi in [-1/np.sqrt(3),0,1/np.sqrt(3)]: ss=s&np.isclose ...` (succeeded in 51ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 12 --ratio 1.25 --rmin 1e-7 --length 16 --output fem_coarse' in /work` (no result line)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 20 --ratio 1.15 --rmin 1e-8 --length 12 --output fem_short' in /work` (succeeded in 2582ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 checks.py' in /work` (succeeded in 52ms:)

## Message 4

The finer mesh changes the phases by less than 0.0003°, and halving the domain length changes them by less than 0.000001°. Direct traction extraction and the energy-integral checks also agree. The combined-load phase is 30.4333°; I’m writing the method and verification details now.
apply patch
patch: completed
/work/method.md
diff --git a/checks.py b/checks.py
new file mode 100644
index 0000000000000000000000000000000000000000..f8b79e6f113ae981a2bed71774ec465a54858228
--- /dev/null
+++ b/checks.py
@@ -0,0 +1,115 @@
+#!/usr/bin/env python3
+"""Compare saved FEM calculations, check them, and write results.json."""
+
+import json
+from pathlib import Path
+
+import numpy as np
+
+from solve import estimate, exact_values
+
+
+def angles(k):
+    return np.degrees(np.angle(np.r_[k, np.sum(k)]))
+
+
+def angle_difference(a, b):
+    return (np.asarray(a) - np.asarray(b) + 180) % 360 - 180
+
+
+def traction_fit(data, material):
+    # Average the two radial Gauss locations after separate extrapolations.
+    estimates = []
+    for xi in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
+        sel = ((data[:, 0].real >= 1e-4) & (data[:, 0].real <= 1e-2)
+               & np.isclose(data[:, 1].real, xi))
+        r = data[sel, 0].real
+        design = np.column_stack([np.ones(len(r)), r, r**2])
+        columns = slice(2, 4) if material == "upper" else slice(4, 6)
+        estimates.append(np.linalg.lstsq(design, data[sel, columns], rcond=None)[0][0])
+    return np.mean(estimates, axis=0)
+
+
+def main():
+    exact = exact_values()
+    runs = {}
+    for name in ("coarse", "medium", "fine", "short"):
+        prefix = Path("fem_" + name)
+        metadata = json.loads(prefix.with_suffix(".json").read_text())
+        with np.load(prefix.with_suffix(".npz")) as data:
+            k = estimate(data["cod"])
+            j = data["j"]
+            selected_j = j[(j[:, 0] > 0.001) & (j[:, 0] < 0.5), 1:]
+            runs[name] = {
+                "nodes": metadata["nodes"],
+                "phases_deg": angles(k).tolist(),
+                "Kabs_cod": np.abs(k).tolist(),
+                "Kabs_relative_errors": (np.abs(k) / exact["Kabs"] - 1).tolist(),
+                "phase_difference_error_deg": float(
+                    angles(k)[0] - angles(k)[1] - exact["phase_difference_abs"]
+                ),
+                "fit_window_change_deg": angle_difference(
+                    angles(estimate(data["cod"], 1e-5, 1e-3)), angles(k)
+                ).tolist(),
+                "J_relative_error_min": (np.min(selected_j, axis=0) / exact["G"][:2] - 1).tolist(),
+                "J_relative_error_max": (np.max(selected_j, axis=0) / exact["G"][:2] - 1).tolist(),
+            }
+            for side in ("upper", "lower"):
+                kt = traction_fit(data["tractions"], side)
+                runs[name][f"traction_{side}_phases_deg"] = angles(kt).tolist()
+                runs[name][f"traction_{side}_Kabs_relative_errors"] = (
+                    np.abs(kt) / exact["Kabs"] - 1
+                ).tolist()
+        assert np.max(np.abs(metadata["force_balance"])) < 1e-12
+        assert np.max(np.abs(metadata["moment_balance"])) < 1e-12
+        assert metadata["residual_max"] < 1e-5
+
+    fine = runs["fine"]
+    assert np.max(np.abs(fine["Kabs_relative_errors"])) < 2e-5
+    assert np.max(np.abs(fine["fit_window_change_deg"])) < 0.01
+    assert abs(fine["phase_difference_error_deg"]) < 0.005
+    for side in ("upper", "lower"):
+        assert np.max(np.abs(angle_difference(fine[f"traction_{side}_phases_deg"], fine["phases_deg"]))) < 0.2
+        assert np.max(np.abs(fine[f"traction_{side}_Kabs_relative_errors"])) < 0.01
+    assert np.max(np.abs(fine["J_relative_error_min"])) < 0.01
+    assert np.max(np.abs(fine["J_relative_error_max"])) < 0.01
+    mesh_change = angle_difference(runs["fine"]["phases_deg"], runs["medium"]["phases_deg"])
+    length_change = angle_difference(runs["short"]["phases_deg"], runs["medium"]["phases_deg"])
+    assert np.max(np.abs(mesh_change)) < 0.01
+    assert np.max(np.abs(length_change)) < 0.01
+
+    # Use the exact magnitudes and fine-mesh phases to form the combined load.
+    psi = fine["phases_deg"][:2]
+    k_exact_mag = np.array(exact["Kabs"]) * np.exp(1j * np.radians(psi))
+    psi_both = float(np.degrees(np.angle(np.sum(k_exact_mag))))
+    results = {
+        "eps": exact["eps"],
+        "caseA": {
+            "G_E2_h_over_P2": exact["G"][0],
+            "Kabs_sqrt_h_over_P": exact["Kabs"][0],
+            "psi_deg": round(psi[0], 4),
+        },
+        "caseB": {
+            "G_E2_h3_over_M2": exact["G"][1],
+            "Kabs_h32_over_M": exact["Kabs"][1],
+            "psi_deg": round(psi[1], 4),
+        },
+        "psi_both_deg": round(psi_both, 4),
+    }
+    report = {
+        "exact": exact,
+        "runs": runs,
+        "fine_minus_medium_phases_deg": mesh_change.tolist(),
+        "length12_minus_length24_phases_deg": length_change.tolist(),
+        "combined_phase_with_exact_magnitudes_deg": psi_both,
+        "all_checks_passed": True,
+    }
+    Path("checks.json").write_text(json.dumps(report, indent=2) + "\n")
+    Path("results.json").write_text(json.dumps(results, indent=2) + "\n")
+    print(json.dumps(report, indent=2))
+    print("Results:")
+    print(json.dumps(results, indent=2))
+
+
+if __name__ == "__main__":
+    main()
diff --git a/method.md b/method.md
new file mode 100644
index 0000000000000000000000000000000000000000..df00216b5f0fd5ed09629cf014ebe53707487af7
--- /dev/null
+++ b/method.md
@@ -0,0 +1,163 @@
+I used exact beam energies for the energy release rates and stress-intensity
+magnitudes, and a two-dimensional plane-strain finite-element solution for the
+phases. All calculations were local, with NumPy and SciPy; no web resources or
+external data were used.
+
+I set h = E2 = 1, so H = 2 and E1 = 9. The two right-hand sides are (P,M) =
+(1,0) and (0,1). In these units the requested phase is simply arg(K), since
+h^(i eps) = 1. This choice also makes the computed G and |K| equal to their
+requested dimensionless normalizations. The crack faces remain traction-free,
+including where this allows interpenetration, as specified in the question.
+
+For the given elastic constants, kappa1 = kappa2 = 5/3, mu1 = 27/8, and
+mu2 = 3/8. Substitution in the prescribed definition gives
+
+\[
+\epsilon=\frac{\ln(2/3)}{2\pi}=-0.0645317762067041.
+\]
+
+Write Ebar_k = E_k/(1-nu_k^2). In a remote separated layer the stress state
+has only sigma11, and its strain energy density is sigma11^2/(2 Ebar_k).
+Because the bonded state far ahead is unstressed, the steady energy release
+rate is the total remote strain energy per unit beam length:
+
+\[
+G=\frac{1}{2\overline E_1}
+  \left(\frac{P^2}{h}+\frac{12M^2}{h^3}\right)
+ +\frac{1}{2\overline E_2}
+  \left(\frac{P^2}{H}
+       +\frac{12[M+P(h+H)/2]^2}{H^3}\right).
+\]
+
+This also follows by evaluating the J integral on a contour far behind and
+ahead of the tip. The rear cross section contributes the beam strain energy;
+the forward cross section and free surfaces contribute zero. Consequently,
+
+\[
+G_A=\frac{287}{162},\qquad G_B=\frac{34}{27},\qquad
+G_{P=M=1}=\frac{815}{162}=G_A+G_B+2.
+\]
+
+For the traction definition of K in the question, the plane-strain
+interface-crack energy relation is
+
+\[
+G=C|K|^2,\qquad
+C=\frac{(\kappa_1+1)/\mu_1+(\kappa_2+1)/\mu_2}
+        {16\cosh^2(\pi\epsilon)}=\frac{64}{135}.
+\]
+
+The reported magnitudes therefore come from the exact expressions
+|K_A| = sqrt(1435/384) = 1.9331267849436742 and
+|K_B| = sqrt(85/32) = 1.629800601300662.
+
+The numerical model in `solve.py` uses nine-node isoparametric quadrilateral
+elements for the displacement formulation. The plane-strain material matrix,
+with engineering shear strain, is
+
+\[
+D_k=E_k\begin{pmatrix}
+1.5&0.75&0\\0.75&1.5&0\\0&0&0.375
+\end{pmatrix}.
+\]
+
+The domain is [-L,L] x [-2,1]. Nodes on y=0 are shared for x>=0 and duplicated
+for x<0. Each left end receives traction (-sigma11,0) from the prescribed
+remote stress distribution. The top, bottom, crack faces, and right end have
+zero prescribed traction. Three scalar displacement constraints at the right
+end remove rigid translations and rotation. The applied forces and moments
+balance, so these constraints introduce no physical loading; numerical
+reactions are included in the reported residual check.
+
+Inside [-1,1] x [-1,1], the mesh consists of geometrically scaled square rings.
+There are `nside` elements along each side of a ring. The first ring connects
+to the common crack-tip node through elements with collapsed inner edges.
+Their Gauss integration points have positive Jacobians. Thus the model has
+an actual crack tip, rather than a small hole. The first positive radius is
+`rmin`; subsequent radii have an approximately constant ratio until radius 1.
+A rectangular mesh covers the remaining domain. The finest calculation used
+nside=28, radius ratio about 1.10, rmin=1e-8, L=24, and 5x5 Gauss integration.
+It had 25,508 elements and 102,685 displacement nodes.
+
+I extracted K primarily from the displacement jump behind the tip. Define
+
+\[
+\delta(r)=[u_2](-r)+i[u_1](-r),\qquad
+S=\frac1{\overline E_1}+\frac1{\overline E_2}.
+\]
+
+The brackets mean upper-layer displacement minus lower-layer displacement.
+The leading crack-opening field, with the traction convention in the task, is
+
+\[
+\delta(r)=\frac{4S}{\cosh(\pi\epsilon)}
+ \frac{K}{1+2i\epsilon}\sqrt{\frac{r}{2\pi}}r^{i\epsilon}.
+\]
+
+I inverted this expression at crack-face nodes, then fitted the resulting
+complex K estimates to K0 + a r + b r^2 over 1e-4 <= r <= 1e-2. The phase of
+K0 is the reported phase for each individual load. The positive real opening
+prefactor affects magnitudes but not phases. As a separate check of the
+complex convention and phase sign, I extracted
+(sigma22+i sigma12) sqrt(2 pi r) r^(-i eps) on the bonded interface, evaluating
+tractions independently from each material.
+
+For M=Ph, linearity gives K_both = K_A + K_B in the chosen units. I used the
+exact magnitudes and the fine-mesh individual phases in this sum. The result
+is psi_both = 30.4332583 degrees, reported as 30.4333 degrees. Adding or
+averaging the individual phase angles would not give this result.
+
+The following checks passed; detailed numbers are saved in `checks.json`.
+
+| Mesh | Nodes | psi A (degrees) | psi B (degrees) | psi of fitted K_A+K_B (degrees) |
+| --- | ---: | ---: | ---: | ---: |
+| nside=12, ratio=1.25, L=16 | 18,257 | 52.250795 | 4.277864 | 30.434086 |
+| nside=20, ratio=1.15, L=24 | 52,749 | 52.249667 | 4.278676 | 30.433515 |
+| nside=28, ratio=1.10, L=24 | 102,685 | 52.249396 | 4.278694 | 30.433307 |
+
+- The medium-to-fine phase changes were below 0.0003 degrees. Reducing L
+  from 24 to 12 on the medium mesh changed all three phases by less than
+  0.000001 degrees.
+- Moving the fine-mesh extrapolation interval to [1e-5,1e-3] changed phases
+  by at most 0.00055 degrees. The alternative [1e-3,0.03] interval is also
+  recorded in the FEM summary files.
+- The fine-mesh opening-derived magnitudes differed from the exact energy
+  values by less than 0.0006 percent.
+- The exact energy cross term requires
+  cos(psi_A-psi_B) = 1/sqrt(G_A G_B), so the magnitude of the phase difference
+  must be 47.97051474 degrees. The fine-mesh phases differ by
+  47.97070200 degrees, an error of 0.00019 degrees.
+- Direct traction extrapolation from the upper material gave phases
+  52.24162 and 4.33705 degrees; from the lower material it gave 52.25725 and
+  4.31626 degrees. These converge toward the opening-derived phases and
+  confirm the prescribed normal-plus-i-shear convention. The largest
+  fine-mesh traction-derived magnitude error was 0.082 percent.
+- J integrals on square contours with radii between 0.001 and 0.5 agreed
+  with the exact energies within 0.087 percent on the fine mesh.
+- Resultant applied forces and moments balanced to approximately 1e-15.
+  The largest fine-mesh linear-system residual, including constrained
+  degrees of freedom, was 1.1e-7 for loads normalized to unity.
+
+To reproduce all saved calculations and regenerate `results.json` and
+`checks.json`, run the following in this folder:
+
+```bash
+OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 12 --ratio 1.25 --rmin 1e-7 --length 16 --output fem_coarse
+OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 20 --ratio 1.15 --rmin 1e-8 --length 24 --output fem_medium
+OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 28 --ratio 1.10 --rmin 1e-8 --length 24 --quadrature 5 --output fem_fine
+OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 solve.py --nside 20 --ratio 1.15 --rmin 1e-8 --length 12 --output fem_short
+OPENBLAS_NUM_THREADS=1 python3 checks.py
+```
+
+The `.npz` files retain crack-opening estimates, interface-traction estimates,
+and contour J integrals, so the checks can be repeated without refactoring
+the stiffness matrices. The individual `fem_*.json` files record parameters,
+fits, force balance, and residuals. The calculation used Python 3.12.3,
+NumPy 1.26.4, and SciPy 1.11.4.
+
+The energies, epsilon, and the energy-derived magnitudes have analytic
+expressions above. The absolute phases have numerical convergence evidence,
+independent traction and energy checks, and domain-size checks; there is no
+rigorous numerical error bound or comparison with an external reference
+phase solution. Effects excluded by the stated linear, traction-free model,
+including crack-face contact, were not evaluated.
diff --git a/solve.py b/solve.py
new file mode 100644
index 0000000000000000000000000000000000000000..caeaf60b53ff5bf981c68a50245993c12a5a7d70
--- /dev/null
+++ b/solve.py
@@ -0,0 +1,315 @@
+#!/usr/bin/env python3
+"""Local, self-contained plane-strain calculation for PUBLIC_TASK.md.
+
+Uses Q2 isoparametric quadrilaterals, geometrically scaled square rings
+around the crack tip, and a conventional mesh away from the tip.  No
+packages other than NumPy and SciPy are needed.
+"""
+
+import argparse
+import json
+from pathlib import Path
+import time
+
+import numpy as np
+from numpy.polynomial.legendre import leggauss
+from scipy.sparse import coo_matrix
+from scipy.sparse.linalg import splu
+
+
+NU = 1.0 / 3.0
+D = np.array([[1.5, 0.75, 0.0], [0.75, 1.5, 0.0], [0.0, 0.0, 0.375]])
+
+
+def exact_values():
+    """Use E2=h=P(or M)=1, H=2, E1=9 throughout."""
+    mu1, mu2 = 9 / (2 * (1 + NU)), 1 / (2 * (1 + NU))
+    kappa = 3 - 4 * NU
+    eps = np.log((kappa / mu1 + 1 / mu2) / (kappa / mu2 + 1 / mu1)) / (2 * np.pi)
+    energy_factor = ((kappa + 1) / mu1 + (kappa + 1) / mu2) / (
+        16 * np.cosh(np.pi * eps) ** 2
+    )
+    ga = (1 - NU**2) / 2 * (1 / 9 + 1 / 2 + 12 * 1.5**2 / 8)
+    gb = (1 - NU**2) / 2 * (12 / 9 + 12 / 8)
+    gab = (1 - NU**2) / 2 * (1 / 9 + 12 / 9 + 1 / 2 + 12 * 2.5**2 / 8)
+    return {
+        "eps": float(eps),
+        "energy_factor": float(energy_factor),
+        "G": [float(ga), float(gb), float(gab)],
+        "Kabs": [float(np.sqrt(ga / energy_factor)), float(np.sqrt(gb / energy_factor))],
+        "phase_difference_abs": float(np.degrees(np.arccos((gab - ga - gb) / (2 * np.sqrt(ga * gb))))),
+    }
+
+
+def shapes(points):
+    """Q2 nodal shapes, local node order (eta=-1,0,1), then xi=-1,0,1."""
+    points = np.atleast_2d(points)
+    x, y = points.T
+    lx = np.array([x * (x - 1) / 2, 1 - x * x, x * (x + 1) / 2]).T
+    ly = np.array([y * (y - 1) / 2, 1 - y * y, y * (y + 1) / 2]).T
+    dx = np.array([x - 0.5, -2 * x, x + 0.5]).T
+    dy = np.array([y - 0.5, -2 * y, y + 0.5]).T
+    shape = np.einsum("qi,qj->qij", ly, lx).reshape(-1, 9)
+    deriv = np.stack(
+        [np.einsum("qi,qj->qij", ly, dx).reshape(-1, 9),
+         np.einsum("qi,qj->qij", dy, lx).reshape(-1, 9)], axis=-1
+    )
+    return shape, deriv
+
+
+class Mesh:
+    def __init__(self, nside=16, ratio=1.22, rmin=1e-7, length=24.0):
+        assert nside % 2 == 0
+        self.coords, self.conn, self.moduli = [], [], []
+        self.lookup, self.rings, self.left_edges = {}, [], []
+
+        # Homothetic square rings.  The innermost elements have a collapsed
+        # edge at the common crack tip; Gaussian quadrature avoids the apex.
+        nr = int(np.ceil(np.log(1 / rmin) / np.log(ratio)))
+        radii = np.r_[0, np.geomspace(rmin, 1.0, nr + 1)]
+        perimeter = []
+        for j in range(nside):
+            perimeter.append((1.0, -1.0 + 2 * j / nside))
+        for j in range(nside):
+            perimeter.append((1.0 - 2 * j / nside, 1.0))
+        for j in range(nside):
+            perimeter.append((-1.0, 1.0 - 2 * j / nside))
+        for j in range(nside):
+            perimeter.append((-1.0 + 2 * j / nside, -1.0))
+        perimeter = np.array(perimeter)
+        for ri, ro in zip(radii[:-1], radii[1:]):
+            ring = []
+            for j, a in enumerate(perimeter):
+                b = perimeter[(j + 1) % len(perimeter)]
+                corners = np.array([ri * a, ro * a, ro * b, ri * b])
+                upper = (a[1] + b[1]) > 0
+                ring.append(self.add_quad(corners, upper))
+            self.rings.append((ri, ro, ring))
+        self.nside = nside
+
+        # Exterior rectangular mesh. Its central square is replaced by rings.
+        step = 2 / nside
+        xfar = [1.0]
+        dx = step
+        while xfar[-1] < length:
+            xfar.append(min(length, xfar[-1] + dx))
+            dx = min(dx * 1.18, 0.8)
+        xs = np.r_[-np.array(xfar[:0:-1]), np.linspace(-1, 1, nside + 1), xfar[1:]]
+        ys = np.linspace(-2, 1, 3 * nside // 2 + 1)
+        for xa, xb in zip(xs[:-1], xs[1:]):
+            for ya, yb in zip(ys[:-1], ys[1:]):
+                if xa >= -1 and xb <= 1 and ya >= -1:
+                    continue
+                corners = np.array([[xa, ya], [xb, ya], [xb, yb], [xa, yb]])
+                ie = self.add_quad(corners, ya >= 0)
+                if xa == -length:
+                    self.left_edges.append(ie)
+        self.coords = np.array(self.coords)
+        self.conn = np.array(self.conn, dtype=np.int32)
+        self.moduli = np.array(self.moduli)
+        self.edofs = (2 * self.conn[:, :, None] + np.arange(2)).reshape(-1, 18)
+        self.length = length
+
+    @staticmethod
+    def key(x, y, upper):
+        # Rounding only serves to identify the same geometric node built
+        # through two adjacent elements. Stored coordinates are unrounded.
+        x, y = round(float(x), 14), round(float(y), 14)
+        crack_face = int(upper) if y == 0.0 and x < 0 else 0
+        return x, y, crack_face
+
+    def node(self, x, y, upper):
+        key = self.key(x, y, upper)
+        if key not in self.lookup:
+            self.lookup[key] = len(self.coords)
+            self.coords.append((float(x), float(y)))
+        return self.lookup[key]
+
+    def add_quad(self, corners, upper):
+        indices = []
+        a, b, c, d = corners
+        for eta in (0, 0.5, 1):
+            for xi in (0, 0.5, 1):
+                pt = (1 - eta) * ((1 - xi) * a + xi * b) + eta * ((1 - xi) * d + xi * c)
+                indices.append(self.node(*pt, upper))
+        ie = len(self.conn)
+        self.conn.append(indices)
+        self.moduli.append(9.0 if upper else 1.0)
+        return ie
+
+
+def geometry(xy, deriv):
+    jac = np.einsum("eni,qnj->eqij", xy, deriv)
+    det = np.linalg.det(jac)
+    if np.any(det <= 0):
+        raise ValueError("Nonpositive element Jacobian")
+    grad = np.einsum("qnj,eqji->eqni", deriv, np.linalg.inv(jac))
+    return grad, det, jac
+
+
+def assemble(mesh, nq=4):
+    ndof = 2 * len(mesh.coords)
+    ne = len(mesh.conn)
+    xg, wg = leggauss(nq)
+    points = np.array([(x, y) for y in xg for x in xg])
+    weights = np.outer(wg, wg).ravel()
+    _, deriv = shapes(points)
+    # Retain all element entries until converting COO to CSR.
+    values = np.empty((ne, 18, 18))
+    for start in range(0, ne, 512):
+        sl = slice(start, min(start + 512, ne))
+        grad, det, _ = geometry(mesh.coords[mesh.conn[sl]], deriv)
+        nx, ny = grad[..., 0], grad[..., 1]
+        weights_e = det * weights * mesh.moduli[sl, None]
+        xx = np.einsum("eqi,eqj,eq->eij", nx, nx, weights_e, optimize=True)
+        yy = np.einsum("eqi,eqj,eq->eij", ny, ny, weights_e, optimize=True)
+        xy = np.einsum("eqi,eqj,eq->eij", nx, ny, weights_e, optimize=True)
+        values[sl, 0::2, 0::2] = 1.5 * xx + 0.375 * yy
+        values[sl, 1::2, 1::2] = 0.375 * xx + 1.5 * yy
+        values[sl, 0::2, 1::2] = 0.75 * xy + 0.375 * xy.transpose(0, 2, 1)
+        values[sl, 1::2, 0::2] = values[sl, 0::2, 1::2].transpose(0, 2, 1)
+    rows = np.broadcast_to(mesh.edofs[:, :, None], values.shape).ravel()
+    cols = np.broadcast_to(mesh.edofs[:, None, :], values.shape).ravel()
+    stiff = coo_matrix((values.ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()
+    force = np.zeros((ndof, 2))
+    # On the left end the outward traction is (-sigma11,0).
+    sg, sd = shapes(np.column_stack([-np.ones(nq), xg]))
+    for ie in mesh.left_edges:
+        xy = mesh.coords[mesh.conn[ie]]
+        pts = sg @ xy
+        jac_edge = np.linalg.norm(sd[:, :, 1] @ xy, axis=1)
+        y = pts[:, 1]
+        if mesh.moduli[ie] == 9:
+            stresses = np.column_stack([-np.ones(nq), -12 * (y - 0.5)])
+        else:
+            stresses = np.column_stack([0.5 + 2.25 * (y + 1), 1.5 * (y + 1)])
+        load = sg.T @ (-stresses * (wg * jac_edge)[:, None])
+        np.add.at(force, mesh.edofs[ie, 0::2], load)
+    return stiff, force
+
+
+def field(mesh, displacement, ie, points):
+    """Returns displacements gradients and stresses, each with both cases."""
+    _, deriv = shapes(points)
+    grad, _, jac = geometry(mesh.coords[mesh.conn[ie]][None], deriv)
+    ue = displacement[mesh.edofs[ie]].reshape(9, 2, 2)
+    # q, displacement component, spatial derivative, load case
+    du = np.einsum("qnj,nic->qijc", grad[0], ue)
+    strain = np.stack([du[:, 0, 0], du[:, 1, 1], du[:, 0, 1] + du[:, 1, 0]], axis=1)
+    stress = mesh.moduli[ie] * np.einsum("ij,qjc->qic", D, strain)
+    return du, strain, stress, jac[0]
+
+
+def extract(mesh, displacement):
+    exact = exact_values()
+    eps = exact["eps"]
+    # Crack opening: [u2]+i[u1] = 4 S / cosh(pi eps) K sqrt(r/(2pi))
+    #                           * r^(i eps)/(1+2i eps), S=sum(1/Ek').
+    compliance_sum = (1 - NU**2) * (1 + 1 / 9)
+    cod_factor = np.sqrt(2 * np.pi) * (1 + 2j * eps) * np.cosh(np.pi * eps) / (4 * compliance_sum)
+    cod_rows, traction_rows, j_rows = [], [], []
+    for ri, ro, ring in mesh.rings[1:]:
+        # Use the ring outer nodes, which belong to both adjoining elements.
+        nt = mesh.lookup[mesh.key(-ro, 0.0, True)]
+        nb = mesh.lookup[mesh.key(-ro, 0.0, False)]
+        delta = displacement[2 * nt + 1] - displacement[2 * nb + 1] + 1j * (
+            displacement[2 * nt] - displacement[2 * nb]
+        )
+        kc = cod_factor * delta / (ro ** (0.5 + 1j * eps))
+        cod_rows.append([ro, *kc])
+        # Tractions on the interface ahead, evaluated from either material.
+        for xi in (-1 / np.sqrt(3), 0.0, 1 / np.sqrt(3)):
+            r = (ri + ro) / 2 + (ro - ri) * xi / 2
+            ku = []
+            for ie, eta in [(ring[mesh.nside // 2], -1.0), (ring[mesh.nside // 2 - 1], 1.0)]:
+                _, _, stress, _ = field(mesh, displacement, ie, [(xi, eta)])
+                k = (stress[0, 1] + 1j * stress[0, 2]) * np.sqrt(2 * np.pi * r) * r ** (-1j * eps)
+                ku.extend(k)
+            traction_rows.append([r, xi, *ku])
+
+        # J integral on the middle square of every ring.
+        if 1e-4 < ro < 0.9:
+            zg, zw = leggauss(5)
+            points = np.column_stack([np.zeros(5), zg])
+            jvalue = np.zeros(2)
+            for ie in ring:
+                du, strain, stress, jac = field(mesh, displacement, ie, points)
+                dx, dy = jac[:, 0, 1], jac[:, 1, 1]
+                energy = 0.5 * np.sum(strain * stress, axis=1)
+                integral = energy * dy[:, None] - (
+                    stress[:, 0] * dy[:, None] - stress[:, 2] * dx[:, None]
+                ) * du[:, 0, 0] - (
+                    stress[:, 2] * dy[:, None] - stress[:, 1] * dx[:, None]
+                ) * du[:, 1, 0]
+                jvalue += zw @ integral
+            j_rows.append([(ri + ro) / 2, *jvalue])
+    return np.array(cod_rows), np.array(traction_rows), np.array(j_rows)
+
+
+def complex_list(z):
+    return [[float(v.real), float(v.imag)] for v in z]
+
+
+def estimate(cod, lower=1e-4, upper=1e-2):
+    sel = (cod[:, 0].real >= lower) & (cod[:, 0].real <= upper)
+    r = cod[sel, 0].real
+    # Leading Williams field plus the next, O(r), correction to K.
+    design = np.column_stack([np.ones(len(r)), r, r**2])
+    return np.linalg.lstsq(design, cod[sel, 1:], rcond=None)[0][0]
+
+
+def run(args):
+    start = time.time()
+    mesh = Mesh(args.nside, args.ratio, args.rmin, args.length)
+    print(f"mesh: {len(mesh.coords)} nodes, {len(mesh.conn)} elements; {time.time()-start:.2f}s", flush=True)
+    stiff, force = assemble(mesh, args.quadrature)
+    n = len(mesh.coords)
+    loads = force.reshape(n, 2, 2)
+    moment = np.einsum("n,nc->c", mesh.coords[:, 0], loads[:, 1]) - np.einsum(
+        "n,nc->c", mesh.coords[:, 1], loads[:, 0]
+    )
+    print(f"force balance: {np.sum(loads, axis=0)}; moment: {moment}", flush=True)
+    p = mesh.lookup[mesh.key(args.length, 0, True)]
+    q = mesh.lookup[mesh.key(args.length, 1, True)]
+    fixed = [2 * p, 2 * p + 1, 2 * q]
+    free = np.setdiff1d(np.arange(2 * n), fixed)
+    print(f"assembled: {time.time()-start:.2f}s, nnz={stiff.nnz}", flush=True)
+    factor = splu(stiff[free][:, free].tocsc())
+    displacement = np.zeros_like(force)
+    displacement[free] = factor.solve(force[free])
+    residual = stiff @ displacement - force
+    print(f"solved: {time.time()-start:.2f}s; max residual {np.max(np.abs(residual)):.3g}", flush=True)
+    cod, tractions, jvals = extract(mesh, displacement)
+    k = estimate(cod)
+    phases = np.degrees(np.angle(np.r_[k, k.sum()]))
+    summary = {
+        "parameters": vars(args),
+        "nodes": n,
+        "elements": len(mesh.conn),
+        "exact": exact_values(),
+        "K_cod": complex_list(k),
+        "Kabs_cod": np.abs(k).tolist(),
+        "phases_deg_cod": phases.tolist(),
+        "residual_max": float(np.max(np.abs(residual))),
+        "force_balance": np.sum(loads, axis=0).tolist(),
+        "moment_balance": moment.tolist(),
+        "elapsed_seconds": time.time() - start,
+    }
+    for lo, hi in [(1e-5, 1e-3), (1e-4, 1e-2), (1e-3, 0.03)]:
+        kk = estimate(cod, lo, hi)
+        summary[f"fit_{lo:g}_{hi:g}"] = {"K": complex_list(kk), "phases": np.degrees(np.angle(np.r_[kk, kk.sum()])).tolist()}
+    prefix = Path(args.output)
+    prefix.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
+    np.savez_compressed(prefix.with_suffix(".npz"), cod=cod, tractions=tractions, j=jvals)
+    print(json.dumps(summary, indent=2), flush=True)
+
+
+if __name__ == "__main__":
+    parser = argparse.ArgumentParser(description=__doc__)
+    parser.add_argument("--nside", type=int, default=16)
+    parser.add_argument("--ratio", type=float, default=1.22)
+    parser.add_argument("--rmin", type=float, default=1e-7)
+    parser.add_argument("--length", type=float, default=24.0)
+    parser.add_argument("--quadrature", type=int, default=4)
+    parser.add_argument("--output", default="fem_default")
+    run(parser.parse_args())

- `/usr/bin/bash -lc "python3 - <<'PY' import json from pathlib import Path r=json.loads(Path('results.json').read_text()) assert set(r)=={'eps','caseA','caseB','psi_both_deg'} assert set(r['caseA'])=={'G_E2_h_over_P2','Kabs_sqrt_h_over_P','psi_deg'} assert set(r['caseB'])=={'G_E2_h3_over_M2','Kabs_h32 ...` (succeeded in 0ms:)
- `/usr/bin/bash -lc 'OPENBLAS_NUM_THREADS=1 python3 checks.py > tmp/verification.log' in /work` (succeeded in 49ms:)

## Message 5

Created [results.json](/work/results.json), [method.md](/work/method.md), and reproducible solver/check code.

Phases: **52.2494°** (A), **4.2787°** (B), **30.4333°** (combined). Energy, mesh-convergence, and domain-size checks passed.

Remaining unverified: rigorous phase-error bounds and comparison with an external reference solution.

