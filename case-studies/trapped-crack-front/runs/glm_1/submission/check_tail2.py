"""Algebraic traction route, corrected Phi_sol assembly:
Phi_m(p) = sum_k basis[m].c[k] * Gamma(3/2+k) (1+ip)^{-3/2-k}  (c includes 1/scale).
"""
import numpy as np
from scipy.special import gamma as G
from m_stiffness import M_closed
import first_order as fo

CAL12 = (1-1j)/np.sqrt(2*np.pi)
ps = np.array([2e3, 4e3, 8e3, 1.6e4, 3.2e4])
Mfit = np.stack([ps**0.5, ps**-0.5, ps**-1.5, ps**-2.5], axis=1)

def Phi_basis(p, bf):
    k = np.arange(len(bf.c))
    vals = np.array([(1+1j*p)**(-1.5-kk) for kk in k])
    return complex((G(1.5+k)*bf.c) @ vals)

for ML in [6, 8]:
    sh = fo.solve_shear(ML=ML, verbose=False)
    basis = sh['basis']; nb = len(basis)
    A12alg = np.zeros((2, 2), complex)
    for jc in range(2):
        s, res = sh['sols'][jc]
        ax, az = s[:nb], s[nb:]
        def Phi_sol(p):
            v = np.zeros(3, complex)
            v[[0, 2][jc]] = np.sqrt(np.pi)*(1+1j*p)**-0.5
            for m in range(nb):
                v[0] += ax[m]*Phi_basis(p, basis[m])
                v[2] += az[m]*Phi_basis(p, basis[m])
            return v
        Tvals = np.array([M_closed(p, 1.0) @ Phi_sol(p) for p in ps])
        for i, comp in enumerate([0, 2]):
            c = np.linalg.lstsq(Mfit, Tvals[:, comp], rcond=None)[0]
            A12alg[i, jc] = CAL12*c[1]
    Ralg = np.zeros((2, 2), complex)
    for i in range(2):
        for j in range(2):
            Ralg[i, j] = fo.S2PI*(fo.DVEC[j]/2.0)*A12alg[i, j]
    print(f"ML={ML}")
    print("  R algebraic (corrected):")
    print(np.array_str(Ralg, precision=8, suppress_small=True))
