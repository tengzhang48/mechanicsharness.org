#!/usr/bin/env python3
"""
Suddenly pressurised crack in a fluid-saturated poroelastic solid.

Model / method (see ANSWER.md):
 * Mode-I crack -L<=x<=L on y=0, plane strain, infinite Biot poroelastic body.
 * Traction-controlled crack: the *total stress* intensity factor is fixed by
   the prescribed crack-face traction (Michel / Bueckner weight-function
   theorem), so it is time independent and modulus independent.
 * Crack volume is obtained exactly from the reciprocity identity
       V^A - V^B = Int  sigma^B_kk(x,y) * alpha(1-2nu) p^A(x,y)/(2G)  dA
   (B = drained elastic reference crack, A = actual poroelastic solution),
   which follows from Betti's theorem for eigenstrain (pore pressure).
   Validated here numerically:  Int (sigma^B_kk)^2 dA = 2 pi L^2 .
"""
import json, numpy as np
from scipy.special import erfc

# ----------------------------- parameters --------------------------------
nu     = 0.2
nu_u   = 0.4
alpha  = 0.89
L      = 1.0

# ----------------------------- formulas ----------------------------------
# Opening traction on the faces:  tau = -sigma_n.  Michel:  K = -sigma_n sqrt(pi L)
K1     = 1.0                      # loading 1: sigma_n=-1  -> K = +sqrt(pi L)
K2     = 0.0                      # loading 2: sigma_n= 0  -> K = 0
# crack volume of an elastic crack with uniform face traction tau: V = pi(1-nu_eff) tau L^2/G
V1_0   = np.pi * (1.0 - nu_u)     # t=0+ : undrained  nu_u
V1_inf = np.pi * (1.0 - nu)       # t->inf: drained  nu
# short-time diffusive growth (1D boundary layer, erfc): V = V1_0 + A sqrt(ct/L^2)
A      = 16.0 * (nu_u - nu) / np.sqrt(np.pi)
# undrained pore pressure next to the faces in loading 1 (sigma_kk=-2 there)
p_faces_0 = 2.0 * (nu_u - nu) / (alpha * (1.0 - 2.0 * nu))
V2_inf = 0.0                      # loading 2: harmonic p -> no stress -> no opening

results = {
    "K1_0": K1,
    "K1_inf": K1,
    "V1_0": V1_0,
    "V1_inf": V1_inf,
    "V1_short_coeff": A,
    "p_faces_0": p_faces_0,
    "K2_0": K2,
    "K2_inf": K2,
    "V2_inf": V2_inf,
    "formulas": {
        "K1_0":   "K = -sigma_n*sqrt(pi*L)  (Michel/weight-function, modulus- & time-"
                  "independent); sigma_n=-H(t) => K1 = +sqrt(pi*L) => K1/sqrt(pi*L) = 1",
        "K1_inf": "same as K1_0:  K1_inf = sqrt(pi*L) => 1  (time independent)",
        "V1_0":   "V = pi(1-nu_u)*(-sigma_n)*L^2/G = pi(1-nu_u) L^2/G  (undrained instantaneous)",
        "V1_inf": "V = pi(1-nu)*L^2/G  (drained long time)",
        "V1_short_coeff":
                  "V(t)=V1_0 + [16(nu_u-nu)/sqrt(pi)] sqrt(ct/L^2); A = 16(nu_u-nu)/sqrt(pi)",
        "p_faces_0":
                  "undrained p from zeta=0: p = -(nu_u-nu)*sigma_kk/[alpha(1-2nu)]; "
                  "sigma_kk=-2 (sigma_xx=sigma_yy=-1 on the faces) => "
                  "p_faces_0 = 2(nu_u-nu)/[alpha(1-2nu)]",
        "K2_0":   "sigma_n=0 on faces => K2 = -sigma_n*sqrt(pi*L) = 0 (and p=1 is a "
                  "measure-zero initial layer, no stress at t=0+)",
        "K2_inf": "at t->inf p is harmonic (p=1 on faces, 0 at infinity) => eigenstrain "
                  "compatible => no stress => K2_inf = 0",
        "V2_inf": "harmonic p => no stress => no opening => V2_inf = 0",
    },
    "reproduce_command": "python3 solve.py",
}

# numeric values for (c)
print("nu=%.3f nu_u=%.3f alpha=%.3f  (K in units sqrt(pi L), V in units L^2/G)" % (nu, nu_u, alpha))
for k in ["K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff",
          "p_faces_0", "K2_0", "K2_inf", "V2_inf"]:
    print(f"  {k:16s} = {results[k]:.10g}")

with open("results.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nwrote results.json")


# ------------------------- numerical validation --------------------------
def sigma_kk(x, y, L=1.0):
    """trace of the stress of the unit-opening-traction elastic crack."""
    z = x + 1j * y
    w = np.sqrt(z * z - L * L)          # principal branch: cut on the crack
    return 2.0 * (np.real(z / w) - 1.0)

def check_identity():
    t, wt = np.polynomial.legendre.leggauss(1200)
    t = 0.5 * (t + 1); mu = 20 * t**2; wmu = wt * 0.5 * (20 * 2 * t)
    q, wq = np.polynomial.legendre.leggauss(1200)
    th = 0.5 * (q + 1) * 2 * np.pi; wth = wq * np.pi
    MU, TH = np.meshgrid(mu, th, indexing="ij")
    D = np.sinh(MU)**2*np.cos(TH)**2 + np.cosh(MU)**2*np.sin(TH)**2
    S = 2.0*(np.cosh(MU)*np.sinh(MU)/D - 1.0)
    jac = np.sinh(MU)**2 + np.sin(TH)**2
    return np.sum(S**2*jac*np.outer(wmu, wth))

val = check_identity()
print("\n[check] Int (sigma_kk)^2 dA = %.6f   (target 2*pi = %.6f)  ratio=%.6f"
      % (val, 2*np.pi, val/(2*np.pi)))

def check_short_time(ct):
    """direct evaluation of V(t)-V^B via the reciprocity integral with a 1D
    erfc boundary-layer pore-pressure profile; returns the implied A."""
    xs = np.linspace(-0.999*L, 0.999*L, 6001)
    h = np.sqrt(ct)
    ys = np.linspace(0.0, 12.0*h, 4001)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    Sk = sigma_kk(X, Y)
    dp = -1.0*erfc(Y/(2.0*h))            # p_und normalised to 1 here
    integ = np.trapezoid(np.trapezoid(Sk*dp, ys, axis=1), xs)      # one side
    integ *= 2.0                            # both faces
    # deltaV/tau ... express the dimensionless character of A
    return integ

for ct in [1e-4, 1e-5, 1e-6]:
    I = check_short_time(ct)
    # deltaV = alpha(1-2nu)/(2G) p_und * I ; p_und = 2(nu_u-nu)/(alpha(1-2nu))
    dV = alpha*(1-2*nu)/(2.0) * (2*(nu_u-nu)/(alpha*(1-2*nu))) * I   # units L^2/G when G=1
    A_num = dV/np.sqrt(ct)
    print("[check] ct=%.0e  A_num=%.6f   A_formula=%.6f" % (ct, A_num, A))