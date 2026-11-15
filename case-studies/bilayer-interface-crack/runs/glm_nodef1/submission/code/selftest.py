"""Self tests for the Q8 implementation: shape functions and a patch test."""
import numpy as np
import sys
sys.path.insert(0, '/work/code')
from fem import q8_shape, assemble, plane_strain_D, elem_stress_at, BilayerCrack

rng = np.random.default_rng(0)

# ---- 1. shape function derivatives vs finite differences
err = 0.0
for _ in range(50):
    xi, eta = rng.uniform(-0.9, 0.9, 2)
    N, dxi, deta = q8_shape(xi, eta)
    h = 1e-6
    Np, _, _ = q8_shape(xi+h, eta); Nm, _, _ = q8_shape(xi-h, eta)
    _, _, dep = q8_shape(xi, eta+h); _, _, dem = q8_shape(xi, eta-h)
    Np2 = q8_shape(xi, eta+h)[0]; Nm2 = q8_shape(xi, eta-h)[0]
    err = max(err, np.abs(dxi-(Np-Nm)/(2*h)).max(),
              np.abs(deta-(Np2-Nm2)/(2*h)).max(),
              abs(N.sum()-1))
    # partition of unity derivative
    err = max(err, abs(dxi.sum()), abs(deta.sum()))
print("shape function max error:", err)
assert err < 1e-8

# ---- 2. patch test: uniform strain in a distorted mesh region
mesh = BilayerCrack(h=1.0, H=2.0, L=6.0, d1=0.2, d2=0.2)
E, nu = 3.0, 0.3
D = plane_strain_D(E, nu)
u = np.zeros(2*len(mesh.coords))
exx, eyy, gxy = 0.003, -0.001, 0.002
u[0::2] = exx*mesh.coords[:, 0] + 0.5*gxy*mesh.coords[:, 1]
u[1::2] = eyy*mesh.coords[:, 1] + 0.5*gxy*mesh.coords[:, 0]
el = mesh.elems[10]
s = elem_stress_at(mesh.coords, el, D, u, *mesh.coords[el[0]])
s_exact = D @ np.array([exx, eyy, gxy])
print("patch test stress:", s, "exact:", s_exact)
assert np.allclose(s, s_exact, atol=1e-10)

# ---- 3. element equilibrium: K*u = f consistency (internal): sum of Ke rows = 0
K = assemble(mesh.coords, mesh.elems, mesh.mats-1, [D, D])
rowsum = np.abs(K @ np.ones(K.shape[0]))
print("max |K*1| (should be ~0):", rowsum.max())
assert rowsum.max() < 1e-9
print("ALL SELF-TESTS PASSED")
