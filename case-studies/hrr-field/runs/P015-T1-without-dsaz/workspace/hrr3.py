"""
HRR plane-strain Mode I, Airy form, solve_bvp with UNKNOWN radial exponent p.

Phi = r^k phi(theta), k = 2 - p.
  sigtilde_rr = A = k phi + phi''
  sigtilde_tt = B = k (k-1) phi
  sigtilde_rt = C = -(k-1) phi'

Compatibility/constitutive ODE identical in structure to hrr2 but with p as a
parameter, so we can check that the solver returns p = 1/(n+1).

BCs (4):
  phi'(0)=0, phi'''(0)=0          (Mode I symmetry)
  A(pi)=0, C(pi)=0                (traction-free faces)
Normalization enforced by the extra condition phi(0)=1 applied via a large
weight? solve_bvp can't take 5 BCs for 4 states + 1 param unless 5 unknowns.
Indeed: 4 state BCs + param p -> S = 4 (+p). To make it square we need 4 BCs
for the 4 states and 1 equation determining p. Instead we use 5 BC components:
  [phi'(0), phi'''(0), phi(0)-1, A(pi), C(pi)]
with unknown vector S = p (size 1).  That is 5 BCs for 4 states + 1 unknown =
5 unknowns -> well posed (BC count = states.size + p.size).
"""
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp

n = 6.5


def rhs(theta, Y, par):
    p = par[0]
    k = 2.0 - p
    phi, d1, d2, d3 = Y
    A = k*phi + d2
    B = k*(k-1.0)*phi
    C = -(k-1.0)*d1
    se = np.sqrt(0.75*(A-B)**2 + 3.0*C**2)
    Ap = k*d1 + d3
    Bp = k*(k-1.0)*d1
    Cp = -(k-1.0)*d2
    sep = (0.75*2.0*(A-B)*(Ap-Bp) + 3.0*2.0*C*Cp)/(2.0*se)
    q = se**(n-1.0)
    qp = (n-1.0)*se**(n-2.0)*sep
    Urp = -(1.0/(2.0*p))*(qp*(A-B) + q*(Ap-Bp))
    Ur = -q*(A-B)/(2.0*p)
    Urpp = (p+1.0)*(q*(B-A)/2.0 - Ur) + 2.0*(qp*C + q*Cp)
    Bpp = k*(k-1.0)*d2
    Cpp = -(k-1.0)*d3
    UmB = A-B
    UmB_p = Ap-Bp
    UmB_pp_noD3 = k*d2 - Bpp
    se2pp = 0.75*(2.0*UmB_p*UmB_p + 2.0*UmB*UmB_pp_noD3) + 3.0*(2.0*Cp*Cp + 2.0*C*Cpp)
    se2pp_d3 = 1.5*UmB
    sepp_noD3 = (se2pp - 2.0*sep*sep)/(2.0*se)
    sepp_d3 = se2pp_d3/(2.0*se)
    qpp_noD3 = (n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep + se**(n-2.0)*sepp_noD3)
    qpp_d3 = (n-1.0)*se**(n-2.0)*sepp_d3
    lhs_noD3 = qpp_noD3*UmB + 2.0*qp*UmB_p + q*(k*d2 - Bpp)
    lhs_d3 = qpp_d3*UmB + q
    d3p = (-2.0*p*Urpp - lhs_noD3)/lhs_d3
    return np.vstack([d1, d2, d3, d3p])


def bc(ya, yb, par):
    p = par[0]
    k = 2.0 - p
    return np.array([
        ya[1],                      # phi'(0)=0
        ya[3],                      # phi'''(0)=0
        ya[0] - 1.0,                # normalization phi(0)=1
        k*yb[0] + yb[2],            # A(pi)=0
        yb[1],                      # C(pi)=0  (C = -(k-1) phi')
    ])


if __name__ == "__main__":
    th = np.linspace(0.0, np.pi, 400)
    y0 = np.zeros((4, th.size))
    y0[0] = 1.0 + 0.1*np.cos(th)
    y0[1] = -0.1*np.sin(th)
    y0[2] = -0.1*np.cos(th)
    p0 = np.array([1.0/(n+1.0)])
    sol = solve_bvp(rhs, bc, th, y0, p=p0, tol=1e-9, max_nodes=40000, verbose=0)
    print("success", sol.success, sol.message, "nodes", sol.x.size)
    if sol.success:
        print("p =", sol.p, " 1/(n+1) =", 1.0/(n+1.0))