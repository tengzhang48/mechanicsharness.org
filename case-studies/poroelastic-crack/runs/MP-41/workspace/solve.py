#!/usr/bin/env python3
"""
MP-41: A suddenly pressurized crack in a fluid-saturated porous solid.

Model / method (see ANSWER.md for the derivation):

The crack is a mode-I (Griffith) crack, half-length L, in an infinite
plane-strain Biot solid described by the constitutive laws given in TASK.md.
Because the material response at t=0+ is UNDRAINED (nu_u) and at t->inf is
DRAINED (nu), the crack opening (mode-I) at those limits is the classical
Griffith result evaluated with the appropriate Poisson ratio, driven by the
total face traction.  The mode-I stress-intensity factor for a crack with
prescribed face traction is independent of the elastic constants, so it is
the same at t=0+ and t->inf.

Elastic (classical) results used (uniform opening pressure sigma, per unit
loading):
    K = sigma * sqrt(pi*L)
    V = pi*(1-nu)*sigma*L^2/G
    (COD) delta(x) = 2*sigma*(1-nu)*sqrt(L^2-x^2)/G

Loadings:
    Loading 1: sigma_n = -H(t) (total normal traction), p = 0 on faces.
    Loading 2: sigma_n = 0 (total normal traction), p = H(t) on faces.
    Loading 1 + Loading 2 = fluid-pressurized crack (sigma_n=-H, p=H).

Short-time growth of V for loading 1:
    The crack faces are drained (p=0) while the interior is undrained; a
    pore-pressure boundary layer of thickness ~sqrt(c t) grows from each
    face.  In this layer the field is 1-D (only y-dependence) and the
    pore pressure obeys  d p/dt = c d2p/dy2  (this is verified below),
    with p(y,0)=p_und and p(0,t)=0, giving p = p_und*erf(y/(2 sqrt(c t))).
    Integrating the vertical-strain change over the layer gives the
    leading-order crack-volume growth.
"""
import math
import numpy as np

# ---- parameters (part c) -------------------------------------------------
nu = 0.2
nu_u = 0.4
alpha = 0.89
L = 1.0
G = 1.0

# ---- quantities ----------------------------------------------------------
# Loading 1
K1_0 = 1.0                      # units of sqrt(pi L)   (= sqrt(pi L))
K1_inf = 1.0                    # units of sqrt(pi L)
V1_0 = math.pi * (1.0 - nu_u)   # units of L^2/G  (undrained)
V1_inf = math.pi * (1.0 - nu)   # units of L^2/G  (drained)

# coefficient A in  V ~ V1_0 + A*sqrt(c t / L^2)
A = 8.0 * (nu_u - nu) / ((1.0 - nu) * math.sqrt(math.pi))

# undrained pore pressure next to a crack face, loading 1 (unit loading)
p_faces_0 = 2.0 * (nu_u - nu) / (alpha * (1.0 - 2.0 * nu))

# Loading 2  (total traction zero on the faces)
K2_0 = 0.0
K2_inf = 0.0
V2_inf = 0.0


# ---- numerical check of the boundary-layer results ----------------------
def boundary_layer_check(c=1.0, p_und=1.0, tmax=0.5, ny=400, ymax=10.0):
    """Independently solve dp/dt = c d2p/dy2 with p(0,t)=0, p(y,0)=p_und,
    p(ymax,t)=p_und, using an implicit stiff ODE integrator (method of lines).
    Return (t, face_value(=0), integral of (p - p_und) dy, analytic value)."""
    from scipy.integrate import solve_ivp
    dy = ymax / ny
    y = np.linspace(0.0, ymax, ny + 1)

    def rhs(t, p):
        d2 = np.zeros_like(p)
        d2[1:-1] = (p[2:] - 2.0 * p[1:-1] + p[:-2]) / dy**2
        d = c * d2
        d[0] = 0.0        # p[0] fixed at 0
        d[-1] = 0.0       # p[-1] fixed at p_und
        return d

    p0 = np.full(ny + 1, p_und)
    p0[0] = 0.0
    sol = solve_ivp(rhs, (0.0, tmax), p0, method="BDF",
                    t_eval=[tmax], rtol=1e-8, atol=1e-10)
    p = sol.y[:, -1]
    p[0] = 0.0
    p[-1] = p_und
    integral = np.trapezoid(p - p_und, y)
    analytic = -p_und * 2.0 * math.sqrt(c * tmax) / math.sqrt(math.pi)
    return tmax, float(p[0]), float(integral), float(analytic)


t_chk, face_chk, int_chk, int_an = boundary_layer_check()

# reproduce the printed values
print("K1_0          =", K1_0)
print("K1_inf        =", K1_inf)
print("V1_0          =", V1_0)
print("V1_inf        =", V1_inf)
print("V1_short_coeff=", A)
print("p_faces_0     =", p_faces_0)
print("K2_0          =", K2_0)
print("K2_inf        =", K2_inf)
print("V2_inf        =", V2_inf)
print("--- checks ---")
print("boundary-layer face value   =", face_chk, "(expect 0)")
print("int(p-p_und)dy  numeric     =", int_chk)
print("int(p-p_und)dy  analytic    =", int_an)
print("relative difference         =", abs(int_chk - int_an) / abs(int_an))

# ---- JSON output ---------------------------------------------------------
import json
results = {
    "K1_0": K1_0,
    "K1_inf": K1_inf,
    "V1_0": V1_0,
    "V1_inf": V1_inf,
    "V1_short_coeff": A,
    "p_faces_0": p_faces_0,
    "K2_0": K2_0,
    "K2_inf": K2_inf,
    "V2_inf": V2_inf,
    "formulas": {
        "K1_0": "sqrt(pi*L)",
        "K1_inf": "sqrt(pi*L)",
        "V1_0": "pi*(1-nu_u)*L^2/G",
        "V1_inf": "pi*(1-nu)*L^2/G",
        "V1_short_coeff": "8*(nu_u-nu)*L^2/((1-nu)*sqrt(pi)*G)",
        "p_faces_0": "2*(nu_u-nu)/(alpha*(1-2*nu))",
        "K2_0": "0",
        "K2_inf": "0",
        "V2_inf": "0",
    },
    "reproduce_command": "python3 solve.py",
}
with open("results.json", "w") as f:
    json.dump(results, f, indent=2)
print("wrote results.json")