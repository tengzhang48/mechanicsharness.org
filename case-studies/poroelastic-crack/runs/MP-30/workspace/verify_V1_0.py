#!/usr/bin/env python3
"""
Independent verification of V1_0 (loading 1, t=0+) via the undrained system.

At t=0+ the fluid has not moved: zeta = 0.  Using the given relations,
zeta=0 gives
    p = -(nu_u - nu)/(alpha(1-2nu)) * sigma_kk                     (3)
Substitute (3) into (1):
    2G e_ij = sigma_ij - nu sigma_kk d_ij
              + alpha(1-2nu) d_ij * [-(nu_u-nu)/(alpha(1-2nu)) sigma_kk]
            = sigma_ij - nu sigma_kk d_ij - (nu_u-nu) sigma_kk d_ij
            = sigma_ij - nu_u sigma_kk d_ij
i.e. Hooke's law with nu_u.  (Already checked symbolically in derive.py.)

Hence the instantaneous crack opening for the uniform face traction -1 is
    w(x,0+) = 4(1-nu_u)/G * P sqrt(L^2-x^2),   P=1
    V1_0    = 2 pi (1-nu_u) L^2 / G.

We verify (i) the algebraic substitution above, and (ii) V1_0 by quadrature of
the opening profile using the undrained nu_u.
"""
import numpy as np
import sympy as sp


def algebraic_check():
    nu, nuu, alpha, G = sp.symbols('nu nu_u alpha G', positive=True)
    s11, s22, p = sp.symbols('s11 s22 p')
    skk = s11 + s22
    # zeta = 0 relation (2):
    p_sol = sp.simplify(-(nuu - nu) / (alpha * (1 - 2 * nu)) * skk)
    lhs = s11 - nu * skk + alpha * (1 - 2 * nu) * p_sol
    return sp.simplify(lhs - (s11 - nuu * skk))


if __name__ == "__main__":
    print("undrained reduction residual (should be 0):", algebraic_check())
    nu, nuu, G, L = 0.2, 0.4, 1.0, 1.0
    x = np.linspace(-L, L, 2000001)
    w = 4 * (1 - nuu) / G * 1.0 * np.sqrt(np.maximum(L**2 - x**2, 0))
    V0num = np.trapezoid(w, x)
    V0ana = 2 * np.pi * (1 - nuu) / G * L**2
    print("V1_0 numeric =", V0num, " analytic =", V0ana,
          " rel err =", abs(V0num - V0ana) / V0ana)