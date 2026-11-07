"""
MP-46 : A suddenly pressurized crack in a fluid-saturated porous solid.

Biot plane-strain poroelasticity (Rice & Cleary 1976 conventions):

    2G e_ij = sigma_ij - nu sigma_kk delta_ij + alpha(1-2nu) delta_ij p
    2G zeta = alpha(1-2nu) sigma_kk + alpha^2 (1-2nu)^2 p /(nu_u - nu)
    sigma_ij,j = 0 ,  q_i = -kappa p,i ,  zeta,t + q_i,i = 0
    c = 2 kappa G (1-nu)(nu_u-nu) / [alpha^2 (1-2nu)^2 (1-nu_u)]

Crack  -L <= x <= L, y = 0.
  Loading 1 : sigma_n = -H(t) , p = 0        (mechanical, drained face)
  Loading 2 : sigma_n =  0    , p = H(t)     (fluid, pressurised face)

This script prints every quantity required by the task, in the units
requested (K in units of sqrt(pi L) ; V in units of L^2/G ; per unit load).

Run:  python3 solve.py
"""
import json
import math

# ----------------------------------------------------------------------
# parameters (part (c))
# ----------------------------------------------------------------------
nu = 0.2
nu_u = 0.4
alpha = 0.89

G = 1.0
L = 1.0

s0 = 1.0    # magnitude of loading 1  (sigma_n = -s0 H(t))
p0 = 1.0    # magnitude of loading 2  (p = p0 H(t))

# ----------------------------------------------------------------------
# Common elastic-crack ingredient (plane strain).
#
#   A 2-D crack |x|<L whose faces carry a uniform *opening* traction
#   P (P>0 opens the crack) has
#          K = P sqrt(pi L) ,      V = pi (1-nu) P L^2 / G .
#   For a crack whose faces are NOT mechanically loaded but carry an
#   effective (skeleton) traction, the same formulae apply with P equal
#   to that effective normal traction.
# ----------------------------------------------------------------------

# ------------------------- Loading 1 ----------------------------------
# sigma_n = -1 , p = 0  ->  opening traction P = +1.
# The total traction on the drained crack faces is constant in time and
# p = 0 there, so the total-stress field (hence K) is time independent;
# only the compliance (and therefore V) relaxes undrained -> drained.
K1_0 = s0 * math.sqrt(math.pi * L)
K1_inf = s0 * math.sqrt(math.pi * L)
V1_0 = math.pi * (1.0 - nu_u) * s0 * L**2 / G      # t = 0+  (undrained)
V1_inf = math.pi * (1.0 - nu) * s0 * L**2 / G      # t -> inf (drained)

# short-time growth: the undrained pore pressure (confined to a boundary
# layer ~ sqrt(c t)) relaxes to 0.  Integrating the erf profile against
# the crack-line weight (1 - x^2/L^2)^{1/2} gives the coefficient below.
V1_short_coeff = 2.0 * math.sqrt(math.pi) * (nu_u - nu) * s0 * L**2 / G

# Undrained pore pressure induced on the crack line by loading 1.
# On the drained-face crack line the elastic stress field is modulus
# independent: sigma_kk = sigma_xx + sigma_yy = -2 (both = -P = -1), so
# at t = 0+   p = -(nu_u-nu) sigma_kk /(alpha(1-2nu)).
p_faces_0 = -1.0 * (nu_u - nu) * (-2.0 * s0) / (alpha * (1.0 - 2.0 * nu))

# ------------------------- Loading 2 ----------------------------------
# sigma_n = 0 , p = 1 on the faces.  The solid *skeleton* feels the
# effective traction  t' = sigma' . n = alpha p n .  For the upper face
# the outward normal is n = (0,-1), so t' = (0, -alpha p): the skeleton
# is pulled *into* the crack  ->  the crack CLOSES.
# The equivalent opening traction is therefore P = -alpha p0, giving
#   K2 = -alpha p0 sqrt(pi L)  (closing, at every time), and
#   V2 = -pi (1-nu) alpha p0 L^2 / G  (drained, t -> inf).
K2_0 = -alpha * p0 * math.sqrt(math.pi * L)
K2_inf = -alpha * p0 * math.sqrt(math.pi * L)
V2_inf = -math.pi * (1.0 - nu) * alpha * p0 * L**2 / G

# ----------------------------------------------------------------------
# report (units: K/sqrt(pi L) , V G/L^2)
# ----------------------------------------------------------------------
RTL = math.sqrt(math.pi * L)
RL2G = L**2 / G

results = {
    "K1_0": K1_0 / RTL,
    "K1_inf": K1_inf / RTL,
    "V1_0": V1_0 / RL2G,
    "V1_inf": V1_inf / RL2G,
    "V1_short_coeff": V1_short_coeff / RL2G,
    "p_faces_0": p_faces_0,
    "K2_0": K2_0 / RTL,
    "K2_inf": K2_inf / RTL,
    "V2_inf": V2_inf / RL2G,
    "formulas": {
        "K1_0": "s0*sqrt(pi*L)",
        "K1_inf": "s0*sqrt(pi*L)",
        "V1_0": "pi*(1-nu_u)*s0*L^2/G",
        "V1_inf": "pi*(1-nu)*s0*L^2/G",
        "V1_short_coeff": "2*sqrt(pi)*(nu_u-nu)*s0*L^2/G",
        "p_faces_0": "2*(nu_u-nu)*s0/(alpha*(1-2*nu))",
        "K2_0": "-alpha*p0*sqrt(pi*L)",
        "K2_inf": "-alpha*p0*sqrt(pi*L)",
        "V2_inf": " -pi*(1-nu)*alpha*p0*L^2/G",
    },
    "reproduce_command": "python3 solve.py",
}

with open("results.json", "w") as fh:
    json.dump(results, fh, indent=2)
    fh.write("\n")

print(json.dumps(results, indent=2))