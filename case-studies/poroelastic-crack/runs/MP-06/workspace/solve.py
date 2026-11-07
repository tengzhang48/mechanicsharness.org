#!/usr/bin/env python3
"""Suddenly pressurized crack in a fluid-saturated poroelastic solid (plane strain).

Analytic results (derived in ANSWER.md), evaluated for nu=0.2, nu_u=0.4, alpha=0.89.
Units: K in units of sigma*sqrt(pi*L); V in units of L^2/G; per unit loading magnitude.
Writes results.json and prints every quantity.
"""
import json
import numpy as np

# ---- material parameters (part c) ----
nu = 0.2          # drained Poisson's ratio
nu_u = 0.4        # undrained Poisson's ratio
alpha = 0.89      # Biot coefficient

gamma = alpha * (1 - 2 * nu) / (1 - nu)   # poroelastic stress-coupling coefficient

res = {}

# ---------------- loading 1: sigma_n = -H(t), p = 0 on crack faces ----------------
# Stress intensity factor is fixed by the traction boundary condition (SIF of a
# uniformly loaded crack is independent of the elastic constants and of time):
res["K1_0"] = 1.0                       # K = 1 * sqrt(pi L)
res["K1_inf"] = 1.0

# Crack volume: undrained (t=0+) and drained (t->inf) elastic plane-strain values
# V = pi*(1-nu_eff)*L^2/G   (Sneddon: V = 2*pi*L^2*sigma*(1-nu^2)/E)
res["V1_0"] = float(np.pi * (1 - nu_u))     # in units L^2/G per unit sigma0
res["V1_inf"] = float(np.pi * (1 - nu))

# Short-time growth: V ~ V1_0 + A*sqrt(c t / L^2)  (diffusion boundary layer at faces)
res["V1_short_coeff"] = float(8 * (nu_u - nu) / (np.sqrt(np.pi) * (1 - nu)))

# Skempton pore pressure induced in the solid adjacent to the faces at t=0+:
# undrained constraint 0 = 2G*zeta with in-plane trace s = sigma_xx+sigma_yy = -2*sigma0
# (pressurized crack carries sigma_xx = sigma_yy = -sigma0 on its faces) =>
# p = -(nu_u-nu)*s/(alpha*(1-2nu)) = 2(nu_u-nu)*sigma0/(alpha*(1-2nu))
res["p_faces_0"] = float(2 * (nu_u - nu) / (alpha * (1 - 2 * nu)))

# ---------------- loading 2: p = H(t) on crack faces, sigma_n = 0 ----------------
# t=0+: no mechanical effect of the face pore pressure (bulk still undrained, zero
# eigenstrain in the bulk, traction-free faces) -> zero stress field:
res["K2_0"] = 0.0

# t->inf: p_ss harmonic (strip potential); the no-crack eigenstress trace on the
# crack plane is sigma_yy^sm = -(gamma/2) p_ss = -gamma/2 (compression: the
# pressurized matrix near the crack is constrained by the far field, like a hot
# inclusion). The traction-free crack relieves this compression by CLOSING, i.e.
# it behaves as a crack with negative internal pressure -gamma/2:
res["K2_inf"] = float(-gamma / 2)                      # K = -(gamma/2) * p0 * sqrt(pi L)
res["V2_inf"] = float(-np.pi * alpha * (1 - 2 * nu) / 2)  # = -pi*(1-nu)*(gamma/2)*L^2/G

# ---------------- general formulas (strings) ----------------
res["formulas"] = {
    "K1": "K1(t) = sigma0*sqrt(pi*L) for all t (SIF fixed by crack-face tractions); "
          "K1_0 = K1_inf = 1 [units sqrt(pi*L) per unit sigma0]",
    "V1_0": "V1(0+) = pi*(1-nu_u)*sigma0*L^2/G   (undrained elastic response, nu->nu_u)",
    "V1_inf": "V1(inf) = pi*(1-nu)*sigma0*L^2/G   (drained elastic response)",
    "V1_short_coeff": "V1 = V1_0 + [8*(nu_u-nu)/(sqrt(pi)*(1-nu))] * (sigma0*L^2/G) "
                      "* sqrt(c*t/L^2) + O(c*t/L^2)",
    "p_faces_0": "p_faces(0+) = 2*(nu_u-nu)*sigma0/(alpha*(1-2*nu))  "
                 "(Skempton value adjacent to the faces; faces carry sigma_xx=sigma_yy=-sigma0)",
    "K2_0": "K2(0+) = 0   (pore pressure on faces has no instantaneous mechanical effect)",
    "K2_inf": "K2(inf) = -(gamma/2)*p0*sqrt(pi*L)  (crack tends to close),  gamma = alpha*(1-2*nu)/(1-nu)",
    "V2_inf": "V2(inf) = -pi*(1-nu)*(gamma/2)*p0*L^2/G = -pi*alpha*(1-2*nu)*p0*L^2/(2*G) (negative = closure/interpenetration in the linear solution)",
    "diffusivity_check": "given c = 2*kappa*G*(1-nu)*(nu_u-nu)/(alpha^2*(1-2nu)^2*(1-nu_u)) "
                         "reproduced by the 2D-trace (indices 1,2) reading of the constitutive laws",
}
res["reproduce_command"] = "python3 solve.py"

if __name__ == "__main__":
    print("nu = %.3f, nu_u = %.3f, alpha = %.3f, gamma = %.6f" % (nu, nu_u, alpha, gamma))
    print("")
    print("--- loading 1 (sigma_n = -H(t), p=0 on faces) ---")
    print("K1_0            = %.6f   [K = K*sigma0*sqrt(pi*L)]" % res["K1_0"])
    print("K1_inf          = %.6f" % res["K1_inf"])
    print("V1_0            = %.6f   [V = V*sigma0*L^2/G]" % res["V1_0"])
    print("V1_inf          = %.6f" % res["V1_inf"])
    print("V1_short_coeff  = %.6f   [V ~ V1_0 + A*sqrt(c t/L^2)]" % res["V1_short_coeff"])
    print("p_faces_0       = %.6f   [p = value*sigma0]" % res["p_faces_0"])
    print("")
    print("--- loading 2 (p = H(t) on faces, sigma_n = 0) ---")
    print("K2_0            = %.6f   [K = K*p0*sqrt(pi*L)]" % res["K2_0"])
    print("K2_inf          = %.6f" % res["K2_inf"])
    print("V2_inf          = %.6f   [V = V*p0*L^2/G]" % res["V2_inf"])
    print("")
    with open("results.json", "w") as f:
        json.dump(res, f, indent=2)
    print("wrote results.json")
