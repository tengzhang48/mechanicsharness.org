"""Steady state (t -> infinity): solve L[Atilde] = 1 on (0,1),
where Atilde(Z) = (E*/2 pi sigma0) dU/dZ is the (negative) climb-dislocation density.
Crack-like singularity built in: Atilde(Z) = R(Z)/sqrt(1-Z), R smooth.
K_inf/(sigma0 sqrt(pi h)) = -pi sqrt(2) * R(1)   (derived below).

Derivation of K relation: crack opening 2u with u = (4K/E*) sqrt(r/(2 pi)),
r = h(1-Z). du/dz = -(2K/E*)/sqrt(2 pi h (1-Z)).
Atilde = (E*/2 pi sigma0) dU/dZ = (E*/2 pi sigma0) du/dz
       = -(K/(pi sigma0 sqrt(2 pi h))) (1-Z)^(-1/2) = R(1)/sqrt(1-Z)
=> K/(sigma0 sqrt(pi h)) = -pi sqrt(2) R(1).
"""
import numpy as np
from kernel import LOperator, chebT


def steady_state(Nc, M=512):
    """Nc unknown Chebyshev coefficients for R."""
    Lop = LOperator(M)
    i = np.arange(1, Nc+1)
    Zc = 0.5*(1.0 - np.cos(np.pi*i/(Nc+0.5)))  # Nc points in (0,1), includes ~1
    Lm = np.zeros((Nc, Nc))
    for k in range(Nc):
        def Hfun(y, k=k):
            return 2.0*chebT(k, 1.0-2.0*np.asarray(y)**2, 0)
        def Hpfun(y, k=k):
            y = np.asarray(y, dtype=float)
            return 2.0*chebT(k, 1.0-2.0*y**2, 1)*(-4.0*y)
        Lm[:, k] = Lop.apply(Zc, Hfun, Hpfun)
    d = np.linalg.solve(Lm, np.ones(Nc))
    R1 = np.sum(d)  # T_k(1)=1
    K_norm = -np.pi*np.sqrt(2.0)*R1
    return K_norm, d, Zc, Lm


if __name__ == "__main__":
    for Nc in [4, 6, 8, 10, 14, 18, 24, 32]:
        K, d, _, _ = steady_state(Nc)
        print(f"Nc={Nc:3d}  R(1)={np.sum(d): .10f}  K_inf/(sigma0 sqrt(pi h)) = {K:.10f}")
    print("classical edge-crack value: 1.1215 (Tada/Wigglesworth)")
