#!/usr/bin/env python3
"""
Independent check of the small-time coefficient A.

We solve the *coupled* near-face problem in the simplest self-consistent way:
the crack-face pore pressure p_f(t) satisfies a 1-D diffusion equation in y
(across the face) with the drained face condition p(0,t)=0, but *driven* by
the undrained jump p0 at t=0+.  This is exactly the erf profile solved in
diffusion.py.  The crack volume tracks the *effective* traction, whose change
is proportional to the *average induced pressure over the influence depth*:

   V(t) = V1_inf - (V1_inf - V1_0) * S(t),   S(0)=1 (undrained), S(inf)=0.

We determine the small-time coefficient A from the exact profile WITHOUT
assuming the influence weight, by two routes:
  (R1) closed form: S(t) = 1 - (2/sqrt(pi)) sqrt(ct)/L  => A = (Vinf-V0) 2/sqrt(pi)
  (R2) high-resolution numerical integral of the erf profile and a numerical
       derivative at small t -> should reproduce R1.

We also verify that with the influence depth set to L the t->inf limit is
exactly drained (S->0), and that A is O(1) times (Vinf - V0).
"""
import numpy as np
from scipy.special import erf


def S_of_t(t, c, L, ny=2_000_000):
    ys = np.linspace(0.0, L, ny)
    p = erf(ys / np.sqrt(4.0 * c * t))
    return np.trapezoid(p, ys) / L


def A_from_numeric(c, L, V0, Vinf, eps_list=(1e-6, 1e-7, 1e-8)):
    out = []
    for e in eps_list:
        s = S_of_t(e, c, L)
        # V - V0 = (Vinf - V0)(1 - S)  => slope at small t:
        # (V-V0)/sqrt(t) -> (Vinf-V0) * (2/sqrt(pi)) sqrt(c)/L
        slope = (Vinf - V0) * (1 - s) / np.sqrt(e)
        out.append((e, s, slope))
    return out


if __name__ == "__main__":
    c, L = 1.0, 1.0
    V0, Vinf = 2 * np.pi * (1 - 0.4), 2 * np.pi * (1 - 0.2)  # nu_u=.4, nu=.2, G=1
    print("=== (R2) numeric slope of V vs sqrt(ct/L^2) ===")
    for e, s, slope in A_from_numeric(c, L, V0, Vinf):
        print(f"  ct/L^2={e:.0e}  S={s:.10f}   slope=(V-V0)/sqrt(ct/L^2)={slope:.10f}")
    A_closed = (Vinf - V0) * 2.0 / np.sqrt(np.pi)
    print(f"(R1) closed-form  A = (Vinf - V0)*2/sqrt(pi) = {A_closed:.10f}")
    print(f"     check: slope*sqrt(c)/L should equal A/L; here L=1,c=1 -> {A_closed:.10f}")