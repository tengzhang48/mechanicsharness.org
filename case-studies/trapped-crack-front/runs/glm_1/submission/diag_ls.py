"""Diagnose the y-sector LS solve: where does the residual live, what does b look like?"""
import numpy as np
from kernels_x import BasisFun, apply_kernels_matrix
from first_order import make_plus_basis, colloc_matrix, B1

ML = 5
basis = make_plus_basis(ML)
Xj = np.logspace(-4, np.log10(25), 40)
Ap = colloc_matrix(basis, Xj, ['yy'])['yy']
A1 = colloc_matrix([B1], Xj, ['yy'])['yy'][:, 0]
w = np.linalg.norm(np.hstack([Ap, A1[:, None]]), axis=1); w[w == 0] = 1
sol = np.linalg.lstsq(Ap/w[:, None], -A1/w, rcond=None)[0]
r = (Ap@sol + A1)
print("coeffs:", sol)
print("cond(Ap weighted):", np.linalg.cond(Ap/w[:, None]))
print("\n  X_j        T[b1]        residual     rel")
for j in range(0, 40, 2):
    print(f"  {Xj[j]:9.5f}  {A1[j]:+12.4e}  {r[j]:+12.4e}  {abs(r[j]/A1[j]):.3f}")
# resulting b profile
b_tot = lambda X: B1.eval(X) + sum(sol[m]*basis[m].eval(X) for m in range(len(basis)))
Xs = np.logspace(-5, 1.5, 20)
print("\n  X        b(X)      b*sqrt(X)   b*X^{3/2}")
for X in Xs:
    print(f"  {X:9.5f}  {b_tot(X):+10.4e}  {b_tot(X)*np.sqrt(X):+10.4e}  {b_tot(X)*X**1.5:+10.4e}")
