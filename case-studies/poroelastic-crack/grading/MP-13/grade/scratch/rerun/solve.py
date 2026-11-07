#!/usr/bin/env python3
"""
MP-13: a suddenly loaded crack in a fluid-saturated poroelastic solid.

Plane-strain Biot solid, crack -L <= x <= L, y = 0.  Two face loadings for t>0:
    Loading 1:  sigma_n = -H(t),  p = 0      (prescribed traction, "dry" crack)
    Loading 2:  sigma_n =  0,     p = +H(t)  (fluid pressure only)

Model / reduction (Rice & Cleary, JGR 81, 409, 1976)
---------------------------------------------------
*  The two 2-D constitutive laws given in the statement reduce to isotropic
    2-D Hooke: the INSTANTANEOUS (undrained) response is governed by (G, nu_u)
    and the long-time (drained) response by (G, nu); the shear modulus G is
    the same in both (checked symbolically in PROCESS.md).
*  A mode-I crack of half-length L under uniform face traction q has
        K  = q sqrt(pi L)
        w  = 2 q (1 - nu_eff)/G  sqrt(L^2 - x^2)
        V  = pi q (1 - nu_eff)/G  L^2
    (the factor is pi, obtained by integrating w).
*  Pore-pressure diffusion makes the effective Poisson ratio relax from nu_u
    (t = 0+) to nu (t -> inf).  For a UNIFORM load the relaxation factor of the
    near-tip field is the same single function of time for every loading,
        Phi(t),   Phi(0+) = (1 - nu_u)/(1 - nu),   Phi(inf) = 1,
    which is exactly the factor Rice & Cleary find for the suddenly pressurised
    crack (per unit fluid pressure).
*  A prescribed TOTAL traction crack (Loading 1) has a time-independent SIF,
    because K = q sqrt(pi L) does not involve the stiffness.  The fluid-pressure
    crack (Loading 2) instead loads the faces through the effective tension
    alpha * p, so its SIF carries the time factor Phi(t).

Reported per unit loading magnitude: K in units of sqrt(pi L), V in L^2/G.
"""

import json
import numpy as np

NU = 0.2
NU_U = 0.4
ALPHA = 0.89
L = 1.0
G = 1.0


def c_coeff(nu, nu_u, alpha, kappa=1.0, G=1.0):
    """The diffusivity exactly as given in the statement (kept for reference)."""
    return 2.0 * kappa * G * (1.0 - nu) * (nu_u - nu) / \
        (alpha ** 2 * (1.0 - 2.0 * nu) ** 2 * (1.0 - nu_u))


def results_core(nu=NU, nu_u=NU_U, alpha=ALPHA, L=L, G=G):
    S = np.sqrt(np.pi * L)
    Phi0 = (1.0 - nu_u) / (1.0 - nu)            # time factor at t = 0+

    # ---------------- Loading 1 : q = -1, p_f = 0 ------------------------
    # prescribed uniform traction -> SIF time independent, = sqrt(pi L)
    K1_0 = S * 1.0
    K1_inf = S * 1.0
    # opening set by the effective Poisson ratio (undrained -> drained)
    V1_0 = np.pi * L ** 2 * (1.0 - nu_u) / G
    V1_inf = np.pi * L ** 2 * (1.0 - nu) / G
    # leading short-time growth of V,  V ~ V1_0 + A sqrt(c t / L^2):
    #   the near-tip boundary layer of width sqrt(c t) relaxes from nu_u to nu.
    A = 8.0 / (3.0 * np.sqrt(np.pi)) * Phi0 * L ** 2 / G

    # pore pressure induced in the solid just outside the faces at t = 0+
    # (the undrained constraint against p = 0 on the faces).  Derived from
    # sigma'_kk = alpha(1 - 2 nu_u)/(nu_u - nu) p  together with the uniform
    # mode-I plane field; general result (see ANSWER.md).
    p_faces_0 = -alpha * (1.0 - 2.0 * nu) * (1.0 - nu) / (2.0 * (2.0 - nu))

    # ---------------- Loading 2 : q = 0, p_f = +1 -------------------------
    # fluid pressure loads the faces through the effective tension alpha * p;
    # drains to the drained value at t -> inf, undrained factor at t = 0+.
    K2_inf = S * alpha * 1.0
    K2_0 = S * alpha * Phi0
    V2_inf = np.pi * alpha * (1.0 - nu) * L ** 2 / G

    return dict(K1_0=K1_0, K1_inf=K1_inf, V1_0=V1_0, V1_inf=V1_inf,
                V1_short_coeff=A, p_faces_0=p_faces_0,
                K2_0=K2_0, K2_inf=K2_inf, V2_inf=V2_inf,
                _Phi0=Phi0)
FORMULAS = {
    "K1_0": "K1_0 = sqrt(pi L)  (prescribed total traction -> time independent)",
    "K1_inf": "K1_inf = sqrt(pi L)",
    "V1_0": "V1_0 = pi L^2 (1 - nu_u)/G  (undrained)",
    "V1_inf": "V1_inf = pi L^2 (1 - nu)/G  (drained)",
    "V1_short_coeff": "A = (8/(3 sqrt(pi))) (1 - nu_u)/(1 - nu) L^2/G, "
                      "V = V1_0 + A sqrt(c t/L^2) + ...  (c as in the statement)",
    "p_faces_0": "p0 = -alpha (1 - 2 nu)(1 - nu) / (2 (2 - nu))",
    "K2_0": "K2_0 = alpha sqrt(pi L) (1 - nu_u)/(1 - nu)",
    "K2_inf": "K2_inf = alpha sqrt(pi L)",
    "V2_inf": "V2_inf = pi L^2 alpha (1 - nu)/G",
}


def main():
    r = results_core()
    out = {
        "K1_0": float(r["K1_0"]),
        "K1_inf": float(r["K1_inf"]),
        "V1_0": float(r["V1_0"]),
        "V1_inf": float(r["V1_inf"]),
        "V1_short_coeff": float(r["V1_short_coeff"]),
        "p_faces_0": float(r["p_faces_0"]),
        "K2_0": float(r["K2_0"]),
        "K2_inf": float(r["K2_inf"]),
        "V2_inf": float(r["V2_inf"]),
        "formulas": FORMULAS,
        "reproduce_command": "python3 solve.py",
    }
    return out


def pretty(out):
    keys = ["K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff",
            "p_faces_0", "K2_0", "K2_inf", "V2_inf"]
    lines = ["=" * 64,
             "MP-13 suddenly loaded poroelastic crack   "
             "nu=%.2f nu_u=%.2f alpha=%.2f" % (NU, NU_U, ALPHA),
             "=" * 64]
    for k in keys:
        lines.append("  %-16s = %.10f" % (k, out[k]))
    lines.append("=" * 64)
    return "\n".join(lines)


if __name__ == "__main__":
    out = main()
    print(pretty(out))
    with open("results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("wrote results.json")