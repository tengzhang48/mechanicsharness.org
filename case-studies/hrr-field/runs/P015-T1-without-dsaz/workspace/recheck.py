"""
Corrected derivation (u ~ r^{p}, eps ~ r^{p-1}).

Airy: Phi = r^k phi, k = 2-p.
  A = sigtilde_rr = k phi + phi''
  B = sigtilde_tt = k(k-1) phi
  C = sigtilde_rt = -(k-1) phi'       (all multiplied by sigma0 K r^{-p})
  se^2 = (3/4)(A-B)^2 + 3 C^2
  s_rr=(A-B)/2, s_tt=(B-A)/2, s_rt=C

Displacement: u_i = (alpha eps0 K^n) r^{p} ubar_i(theta)   [since u~r^p]
  eps_rr     = p r^{p-1} Urbar
  eps_tt     = r^{p-1}(Urbar + Utbar')
  2 eps_rt   = r^{p-1}(Urbar' - p Utbar - Utbar)
Constitutive: eps = (alpha eps0 K^n) r^{p-1} q s   with q = se^{n-1}
  =>  p Urbar              = q (A-B)/2
      Urbar + Utbar'       = q (B-A)/2
      Urbar' - p Utbar - Utbar = 2 q C

From eq1: Urbar = q(A-B)/(2p).
From eq3: (1-p)? : Urbar' - p Utbar - Utbar = 2 q C
   => Urbar' - (p+1) Utbar = 2 q C  => Utbar = (Urbar' - 2 q C)/(p+1)
From eq2: Utbar' = q(B-A)/2 - Urbar
Compatibility: d/dtheta(Utbar) must equal eq2.
"""
import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from scipy.optimize import fsolve

n = 6.5
p = 1.0/(n+1.0)
k = 2.0 - p


def rhs(theta, Y):
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
    # Urbar = q(A-B)/(2p)
    Ur = q*(A-B)/(2.0*p)
    # Urbar' = (1/2p)(qp(A-B)+q(Ap-Bp))
    Urp = (1.0/(2.0*p))*(qp*(A-B) + q*(Ap-Bp))
    # Utbar = (Urp - 2 q C)/(p+1)
    Ut = (Urp - 2.0*q*C)/(p+1.0)
    # eq2: Ut' = q(B-A)/2 - Ur  ; differentiate Ut:
    # Ut' = (Urp' - 2(qp C + q Cp))/(p+1)
    # => Urp' = (p+1)(q(B-A)/2 - Ur) + 2(qp C + q Cp)
    Urpp = (p+1.0)*(q*(B-A)/2.0 - Ur) + 2.0*(qp*C + q*Cp)
    # Urp' = (1/2p)[ qpp(A-B) + 2 qp(Ap-Bp) + q(App-Bpp) ]
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
    # (1/2p)[ qpp UmB + 2 qp UmB_p + q(k*d2 + d3p - Bpp) ] = Urpp
    lhs_noD3 = qpp_noD3*UmB + 2.0*qp*UmB_p + q*(k*d2 - Bpp)
    lhs_d3 = qpp_d3*UmB + q
    # (1/2p)(lhs_noD3 + lhs_d3*d3p) = Urpp
    d3p = (2.0*p*Urpp - lhs_noD3)/lhs_d3
    return [d1, d2, d3, d3p]


def endvals(nn, u):
    global n, p, k
    n = nn; p = 1.0/(nn+1.0); k = 2.0 - p
    Y0 = [1.0, 0.0, u[0], u[1]]
    s = solve_ivp(rhs, [0.0, np.pi], Y0, rtol=1e-8, atol=1e-10,
                  method='Radau', max_step=0.03)
    if not s.success:
        return np.array([1e6, 1e6])
    phi, d1, d2, d3 = s.y[:, -1]
    return np.array([k*phi + d2, -(k-1.0)*d1])


if __name__ == "__main__":
    # solve at n=6.5 directly with fsolve, seeded
    for seed in [[-1.1, 0.1], [-1.0, 0.2], [0.7, 0.0]]:
        sol, info, ier, msg = fsolve(lambda x: endvals(6.5, x), seed,
                                     xtol=1e-11, full_output=True)
        r = endvals(6.5, sol)
        print("seed", seed, "phi2=%.6f phi3=%.6f A(pi)=%.2e C(pi)=%.2e ier=%d"
              % (sol[0], sol[1], r[0], r[1], ier))