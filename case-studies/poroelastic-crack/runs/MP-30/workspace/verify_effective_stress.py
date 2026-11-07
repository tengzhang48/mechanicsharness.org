#!/usr/bin/env python3
"""
Identify the effective stress form of the *given* constitutive law, and use it
to get V2_inf rigorously (loading 2).

Given (2-D indices i,j in {1,2}; sigma_kk = sigma_11+sigma_22):
  2G e_ij = sigma_ij - nu sigma_kk d_ij + alpha (1-2nu) d_ij p     (1)
  2G zeta = alpha(1-2nu) sigma_kk + alpha^2 (1-2nu)^2 p/(nu_u-nu)  (2)

Seek sigma^eff_ij = sigma_ij + beta p d_ij (2-D). Then
  sigma^eff_kk = sigma_kk + 2 beta p
  sigma^eff_ij - nu sigma^eff_kk d_ij
     = sigma_ij - nu sigma_kk d_ij + beta(1-2nu) p d_ij
Match with (1): beta(1-2nu) = alpha(1-2nu)  =>  beta = alpha.

=>  sigma^eff_ij = sigma_ij + alpha p d_ij   (standard Biot effective stress)
    and (1) is Hooke's law in sigma^eff with the *drained* nu.

A crack with faces carrying total traction sigma_n and pore pressure p opens
as an elastic crack (drained nu, G) under uniform effective face tension
    sigma^eff_n = sigma_n + alpha p.

Loading 1 faces: sigma_n=-1, p=0 -> sigma^eff_n=-1
    => V1_inf = 2 pi (1-nu)/G * L^2
Loading 2 faces: sigma_n=0,  p=1 -> sigma^eff_n= alpha
    => V2_inf = 2 pi (1-nu)/G * alpha * L^2 = alpha V1_inf
"""
import numpy as np
import sympy as sp


def verify_effective_stress_2d():
    G, nu, alpha, p = sp.symbols('G nu alpha p', positive=True)
    s11, s22 = sp.symbols('s11 s22')
    b = alpha
    e11 = s11 + b * p
    e22 = s22 + b * p
    ekk = e11 + e22                      # 2-D trace of sigma^eff
    lhs_11 = e11 - nu * ekk
    skk = s11 + s22                      # 2-D trace of sigma
    rhs_11 = s11 - nu * skk + alpha * (1 - 2 * nu) * p
    return sp.simplify(lhs_11 - rhs_11)


if __name__ == "__main__":
    print("2-D effective-stress residual (should be 0):",
          verify_effective_stress_2d())
    nu, nuu, alpha, G, L = 0.2, 0.4, 0.89, 1.0, 1.0
    V1inf = 2 * np.pi * (1 - nu) / G * L * L
    V2inf = 2 * np.pi * (1 - nu) / G * alpha * L * L
    print("V1_inf =", V1inf)
    print("V2_inf =", V2inf, " (= alpha*V1_inf =", alpha * V1inf, ")")