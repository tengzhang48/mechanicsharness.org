import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# =====================================================================
#  Plane-strain power-law (Ramberg-Osgood) V-notch leading asymptote.
#
#  sigma_ij = r^{-s} Sigma_ij(theta),  Airy Phi = r^{2-s} f(theta):
#     Sigma_rr = f'' + (2-s) f
#     Sigma_tt = (2-s)(1-s) f
#     Sigma_rt = (s-1) f'
#  Incompressible plastic:  eps_ij = (3/2) B sigma_e^{n-1} s_ij.
#
#  Reduced nonlinear eigenvalue ODE (C>0 overall factor, cancels -> C=1):
#     D = f'' + s(2-s) f
#     R = (s-1) f'
#     g = sqrt( (3/4) D^2 + 3 R^2 )               # = sigma_e r^s /(unit)
#     P = g^{n-1} D
#     ODE:   P'' + s n (2-s) P - 4 (g^{n-1} R)' = 0
#     V (u_theta ampl.):  V' = -(2-s n)/(1-s n) P
#
#  BCs (mode I, symmetric about theta=0, traction-free at +-alpha):
#     theta=0 : f'(0)=0 , V(0)=0
#     theta=alpha : Sigma_tt=0 , Sigma_rt=0
#  Symmetry also makes Sigma_rt(-)= ... using theta in [0,alpha] with f'(0)=0.
# =====================================================================

def _P(jet, s, n):
    f0, f1, f2, f3, f4 = jet
    D = f2 + s*(2.0 - s)*f0
    R = (s - 1.0)*f1
    g = np.sqrt(0.75*D*D + 3.0*R*R)
    return g**(n - 1.0)*D

def _GR(jet, s, n):
    f0, f1, f2, f3, f4 = jet
    D = f2 + s*(2.0 - s)*f0
    R = (s - 1.0)*f1
    g = np.sqrt(0.75*D*D + 3.0*R*R)
    return g**(n - 1.0)*R

def dP_dth(jet, s, n):
    # (d/dth) P with df0=f1, df1=f2, df2=f3, df3=f4 (P has no f4 dependence)
    f0, f1, f2, f3, f4 = jet
    # complex-step partials
    h = 1e-30
    def part(func, k):
        j1 = list(jet); j1[k] = j1[k] + 1j*h
        return np.imag(func(j1, s, n))/h
    return (part(_P,0)*f1 + part(_P,1)*f2 + part(_P,2)*f3)

def dP2_dth2(jet, s, n):
    # P'' : differentiate P' w.r.t theta; P' depends on f0..f3 with
    # df0=f1..df3=f4.  Use complex-step partials of P' itself.
    f0, f1, f2, f3, f4 = jet
    h = 1e-30
    def Pprime(j):
        return dP_dth(j, s, n)
    def part(k):
        j1 = list(jet); j1[k] = j1[k] + 1j*h
        return np.imag(Pprime(j1))/h
    # d/dth P' = sum_k dP'/df_k * df_k/dth ; df0=f1,df1=f2,df2=f3,df3=f4
    return part(0)*f1 + part(1)*f2 + part(2)*f3 + part(3)*f4

def dGR_dth(jet, s, n):
    f0, f1, f2, f3, f4 = jet
    h = 1e-30
    def part(k):
        j1 = list(jet); j1[k] = j1[k] + 1j*h
        return np.imag(_GR(j1, s, n))/h
    return part(0)*f1 + part(1)*f2 + part(2)*f3

def make_rhs(n, s):
    def rhs(th, y):
        f0, f1, f2, f3, V = y
        # ODE: P'' + s n (2-s) P - 4 (g^{n-1} R)' = 0, linear in f4??
        # P'' depends on f4 through part(3)*f4 term. Solve for f4:
        # P'' = A + B*f4  where A = part(0)*f1+part(1)*f2+part(2)*f3,
        # B = part(3).  Then f4 = (4 (g^{n-1} R)' - s n (2-s) P - A)/B.
        jet = [f0, f1, f2, f3, 0.0]
        h = 1e-30
        def Pprime(j):
            return dP_dth(j, s, n)
        def partPp(k):
            j1 = list(jet); j1[k] = j1[k] + 1j*h
            return np.imag(Pprime(j1))/h
        A = partPp(0)*f1 + partPp(1)*f2 + partPp(2)*f3
        B = partPp(3)
        Pv = _P(jet, s, n)
        GRp = dGR_dth(jet, s, n)
        f4 = (4.0*GRp - s*n*(2.0 - s)*Pv - A)/B
        Vp = -(2.0 - s*n)/(1.0 - s*n)*Pv
        return [f1, f2, f3, f4, Vp]
    return rhs

def integrate(s, n, alpha, f2_0, f3_0, f0_0=1.0):
    y0 = [f0_0, 0.0, f2_0, f3_0, 0.0]
    rhs = make_rhs(n, s)
    return solve_ivp(rhs, [0.0, alpha], y0, rtol=1e-10, atol=1e-12,
                     method='RK45', dense_output=True, max_step=0.02)

def bcs(x, n, alpha):
    s, f2_0, f3_0 = x
    sol = integrate(s, n, alpha, f2_0, f3_0)
    f0, f1, f2, f3, V = sol.y[:, -1]
    Tt = (2.0 - s)*(1.0 - s)*f0
    Rt = (s - 1.0)*f1
    return np.array([Tt, Rt, V])

def solve(n, alpha, s0, f2_0=0.5, f3_0=0.0):
    x, info, ier, msg = fsolve(bcs, [s0, f2_0, f3_0], args=(n, alpha),
                               full_output=True, xtol=1e-12)
    return x, bcs(x, n, alpha), ier, msg

if __name__ == "__main__":
    alpha = 2.0*np.pi/3.0
    print("== crack limit sanity: alpha->pi should give s=1/(n+1) ==")
    for nn in [4.0]:
        for al in [np.pi*0.999, alpha]:
            found = []
            for s0 in np.linspace(0.02, 0.9, 45):
                for f2g in [-2,-1,0,1,2]:
                    try:
                        x, r, ier, msg = solve(nn, al, s0, f2_0=f2g)
                        if np.linalg.norm(r) < 1e-8 and x[0] > 1e-4:
                            found.append(round(x[0],6))
                    except Exception:
                        pass
            found = sorted(set(found))
            print(f" n={nn} alpha={al:.4f}: candidate s = {found}  (target crack {1/(nn+1):.4f})")