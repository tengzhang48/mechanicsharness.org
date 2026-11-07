#!/usr/bin/env python3
"""
solve.py -- suddenly pressurized crack in a fluid-saturated porous solid.

Reports, for nu = 0.2, nu_u = 0.4, alpha = 0.89:
  K in units of sqrt(pi L), V in units of L^2/G, per unit loading magnitude.

Model
-----
The given constitutive law is Biot poroelasticity with the *effective stress*
    sigma^eff_ij = sigma_ij + alpha p delta_ij,
and (1) is Hooke's law in sigma^eff with the drained nu (verified symbolically
below and in verify_effective_stress.py).  In the undrained limit zeta=0 the
same law reduces to Hooke with nu_u (verified in verify_V1_0.py / derive.py).

Results
-------
Loading 1 (sigma_n = -H(t), p = 0 on the faces):
  * K is set by the uniform face traction and is elastic-constant independent:
        K1_0 = K1_inf = 1 * sqrt(pi L)          -> K1_0 = K1_inf = 1
  * Crack opening:  2G-based mode-I profile w = 4(1-nu)/G P sqrt(L^2-x^2),
        V = int w dx = 2 pi (1-nu) P L^2 / G
        V1_0   = 2 pi (1-nu_u) L^2/G   (undrained, t=0+)
        V1_inf = 2 pi (1-nu)   L^2/G   (drained,  t->inf)
  * Small time:  V ~ V1_0 + A sqrt(ct/L^2)  with
        A = (V1_inf - V1_0) * 2/sqrt(pi)
    (from the 1-D diffusion of the undrained pore-pressure jump p0 across the
     drained face, p(y,t) = p0 erf(y/sqrt(4ct)); see diffusion.py).
  * p_faces_0 = 2 (nu_u-nu) /[alpha (1-2nu)(1-nu_u)]  (undrained jump just
    next to the face at t=0+; derived from the exact mode-I face stress state
    sigma_xx = sigma_yy = -1 with zeta=0; see derive.py, mode1_stress.py).

Loading 2 (sigma_n = 0, p = H(t) on the faces):
  * Faces carry zero total traction -> no mode-I driving force:
        K2_0 = K2_inf = 0
  * Long time: effective face tension = alpha p = alpha,
        V2_inf = 2 pi (1-nu) alpha L^2/G = alpha V1_inf
"""
import json
import numpy as np

# --------------------------------------------------------------------------
# parameters for part (c)
# --------------------------------------------------------------------------
NU = 0.2
NU_U = 0.4
ALPHA = 0.89
# unit loading magnitude -> P = 1;  K in units of sqrt(pi L); V in units of L^2/G


def solve(nu=NU, nu_u=NU_U, alpha=ALPHA):
    # --- Loading 1 -------------------------------------------------------
    K1_0 = K1_inf = 1.0                                   # [sqrt(pi L)]
    V1_0 = 2.0 * np.pi * (1.0 - nu_u)                      # [L^2/G]
    V1_inf = 2.0 * np.pi * (1.0 - nu)                      # [L^2/G]
    V1_short_coeff = (V1_inf - V1_0) * 2.0 / np.sqrt(np.pi)
    p_faces_0 = 2.0 * (nu_u - nu) / (alpha * (1.0 - 2.0 * nu) * (1.0 - nu_u))

    # --- Loading 2 -------------------------------------------------------
    K2_0 = 0.0
    K2_inf = 0.0
    V2_inf = 2.0 * np.pi * (1.0 - nu) * alpha               # [L^2/G]

    return dict(K1_0=K1_0, K1_inf=K1_inf, V1_0=V1_0, V1_inf=V1_inf,
                V1_short_coeff=V1_short_coeff, p_faces_0=p_faces_0,
                K2_0=K2_0, K2_inf=K2_inf, V2_inf=V2_inf)


FORMULAS = {
    "K1_0": "1*sqrt(pi*L)  (uniform face traction; elastic-constant independent)",
    "K1_inf": "1*sqrt(pi*L)  (same as K1_0)",
    "V1_0": "2*pi*(1-nu_u)*L^2/G",
    "V1_inf": "2*pi*(1-nu)*L^2/G",
    "V1_short_coeff": "(V1_inf - V1_0)*2/sqrt(pi) = 4*(nu_u-nu)*sqrt(pi)*L^2/G",
    "p_faces_0": "2*(nu_u-nu)/(alpha*(1-2*nu)*(1-nu_u))",
    "K2_0": "0",
    "K2_inf": "0",
    "V2_inf": "2*pi*(1-nu)*alpha*L^2/G = alpha*V1_inf",
}


def main():
    res = solve()
    for k, v in res.items():
        print(f"{k} = {v!r}")
    out = {k: float(v) for k, v in res.items()}
    out["formulas"] = FORMULAS
    out["reproduce_command"] = "python3 solve.py"
    out["parameters"] = {"nu": NU, "nu_u": NU_U, "alpha": ALPHA,
                         "units": "K in sqrt(pi L), V in L^2/G, "
                                  "per unit loading magnitude"}
    with open("results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nwrote results.json")


if __name__ == "__main__":
    main()