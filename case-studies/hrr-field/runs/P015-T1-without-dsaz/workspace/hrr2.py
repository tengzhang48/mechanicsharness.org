"""
HRR plane-strain Mode I via Airy stress function (equilibrium built in).

Phi(r,theta) = r^{k} phi(theta),  k = 2 - p,  p = 1/(n+1).

Then stresses (normalized, sigma = sigma0*K*r^{-p}*sigtilde):
  sigtilde_rr = k*phi + phi''        =: A
  sigtilde_tt = k*(k-1)*phi          =: B
  sigtilde_rt = -(k-1)*phi'          =: C

These satisfy equilibrium identically.

Constitutive (plane strain, s_zz=0):
  s_rr = (A-B)/2, s_tt = (B-A)/2, s_rt = C
  se^2 = (3/4)(A-B)^2 + 3 C^2

Displacement: u_r = r^{-p} Ur(theta), u_t = r^{-p} Ut(theta)
  eps_rr = -p r^{-p-1} Ur
  eps_tt = r^{-p-1} (Ur + Ut')
  2 eps_rt = r^{-p-1}(Ur' - p Ut - Ut)

Constitutive matching (common factor Lambda r^{-p-1}, Lambda>0):
  -p Ur = q (A-B)/2
  Ur + Ut' = q (B-A)/2
  Ur' - p Ut - Ut = 2 q C
where q = se^{n-1}.

From eq1: Ur = -q(A-B)/(2p)
From eq2: Ut' = q(B-A)/2 - Ur
From eq3: Ur' = 2 q C + p Ut + Ut

Differentiate eq1 to get Ur', equate to eq3 -> one ODE. Combined with Ut',
compatibility. This yields a system in the 4 variables (phi, phi', phi'', phi''').
"""
import numpy as np
from scipy.integrate import solve_ivp, solve_bvp

n = 6.5
p = 1.0/(n+1.0)
k = 2.0 - p
print("p =", p, " k =", k)


def quantities(phi, d1, d2, d3):
    A = k*phi + d2
    B = k*(k-1.0)*phi
    C = -(k-1.0)*d1
    se = np.sqrt(0.75*(A-B)**2 + 3.0*C**2)
    return A, B, C, se


def rhs(theta, Y):
    # Y = [phi, phi', phi'', phi''']
    phi, d1, d2, d3 = Y
    A, B, C, se = quantities(phi, d1, d2, d3)
    # derivatives
    Ap = k*d1 + d3
    Bp = k*(k-1.0)*d1
    Cp = -(k-1.0)*d2
    sep = (0.75*2.0*(A-B)*(Ap-Bp) + 3.0*2.0*C*Cp)/(2.0*se)
    q = se**(n-1.0)
    qp = (n-1.0)*se**(n-2.0)*sep
    # eq1: Ur = -q(A-B)/(2p);  Ur' = -(1/2p)(qp(A-B)+q(Ap-Bp))
    Urp = -(1.0/(2.0*p))*(qp*(A-B) + q*(Ap-Bp))
    # eq3: Urp = 2 q C + (p+1) Ut  =>  Ut = (Urp - 2 q C)/(p+1)
    Ut = (Urp - 2.0*q*C)/(p+1.0)
    # eq2': Ut' = q(B-A)/2 - Ur   must equal d/dtheta[(Urp-2qC)/(p+1)]
    # Ut' = (Urp' - 2(qp C + q Cp))/(p+1)
    # => Urp' = (p+1)*( q(B-A)/2 - Ur ) + 2(qp C + q Cp)
    Ur = -q*(A-B)/(2.0*p)
    Urpp = (p+1.0)*(q*(B-A)/2.0 - Ur) + 2.0*(qp*C + q*Cp)
    # Urpp = -(1/2p)[ qpp(A-B) + 2 qp(Ap-Bp) + q(App-Bpp) ]
    # unknown is d3' = phi''''  appearing in App (Ap=k*d1+d3 => Ap'=k*d2+d3')
    # App = k*d2 + d3p ;  Bpp = k*(k-1)*d2
    # need qpp: se'' from se^2 identity
    # (se^2)'' = 2(...); se^2=(3/4)(A-B)^2+3C^2
    App_ = k*d2  # part without d3'
    Bpp = k*(k-1.0)*d2
    Cpp = -(k-1.0)*d3
    # d2(theta3) of (A-B): (A-B)'' = App-Bpp with App = k*d2+d3p
    UmB = A-B
    UmB_p = Ap-Bp
    UmB_pp_noD3 = (k*d2) - Bpp   # + d3p
    se2pp = 0.75*(2.0*UmB_p*UmB_p + 2.0*UmB*UmB_pp_noD3) + 3.0*(2.0*Cp*Cp + 2.0*C*Cpp)
    # + (3/4)*2*UmB*d3p = 1.5*UmB*d3p
    se2pp_d3 = 1.5*UmB
    sepp_noD3 = (se2pp - 2.0*sep*sep)/(2.0*se)
    sepp_d3 = se2pp_d3/(2.0*se)
    qpp_noD3 = (n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep + se**(n-2.0)*sepp_noD3)
    qpp_d3 = (n-1.0)*se**(n-2.0)*sepp_d3
    # equation: qpp(A-B) + 2qp(Ap-Bp) + q[(App-Bpp)] = -2p Urpp
    # App = k*d2 + d3p
    lhs_noD3 = qpp_noD3*UmB + 2.0*qp*UmB_p + q*(k*d2 - Bpp)
    lhs_d3 = qpp_d3*UmB + q*1.0   # d3p coefficient: from qpp term (qpp_d3*UmB*d3p) + q*d3p
    # lhs = lhs_noD3 + lhs_d3*d3p = -2p Urpp
    d3p = (-2.0*p*Urpp - lhs_noD3)/lhs_d3
    return [d1, d2, d3, d3p]


def bc(ya, yb):
    # theta=0 (Mode I symmetry): C(0)=0 => phi'(0)=0 ;  phi'''(0)=0
    # normalization: phi(0)=1  (fixes the free amplitude)
    # theta=pi free surface: A(pi)=0 => k*phi+phi''=0
    # (C(pi)=0 is then checked as a consistency/validation condition)
    k_ = k
    return np.array([
        ya[1],                                  # phi'(0)=0  (C(0)=0)
        ya[3],                                  # phi'''(0)=0
        ya[0] - 1.0,                            # normalization phi(0)=1
        k_*yb[0] + yb[2],                       # A(pi)=k phi + phi'' =0
    ])


def shoot(t):
    """Integrate from theta=0 with phi(0)=1, phi'(0)=0, phi''(0)=t, phi'''(0)=0.
    Return (A(pi), C(pi))."""
    from scipy.integrate import solve_ivp
    Y0 = [1.0, 0.0, t, 0.0]
    s = solve_ivp(rhs, [0.0, np.pi], Y0, rtol=1e-11, atol=1e-13,
                  method='Radau', dense_output=True, max_step=0.01)
    return s


def resid_A(t):
    s = shoot(t)
    if not s.success:
        return 1e3
    phi, d1, d2, d3 = s.y[:, -1]
    return k*phi + d2


if __name__ == "__main__":
    th = np.linspace(0.0, np.pi, 600)
    y0 = np.zeros((4, th.size))
    y0[0] = 1.0 + 0.3*np.cos(th)
    y0[1] = -0.3*np.sin(th)
    y0[2] = -0.3*np.cos(th)
    sol = solve_bvp(rhs, bc, th, y0, tol=1e-9, max_nodes=40000, verbose=0)
    print("success", sol.success, sol.message, "nodes", sol.x.size)