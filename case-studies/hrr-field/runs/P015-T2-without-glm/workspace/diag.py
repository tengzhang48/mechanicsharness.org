import numpy as np
from scipy.integrate import solve_ivp
import solve as S

# Test: does rhs_maker's f'''' agree with the ODE (GD)'' + a GD - 4 (GT)' = 0
# solved for f'''' (via numeric differentiation of the composite terms)?
n, p, alpha = 4.0, 0.15, 2*np.pi/3
c0 = p*(2-p); a = n*p*(2-n*p)/(1-n*p)

def ftest(th):   # smooth test function
    return 1.0 + 0.3*np.cos(1.7*th) - 0.2*np.cos(0.9*th)**2

h = 1e-3
th0 = 0.7
def derivs(th):
    # 4th-order central differences
    f  = ftest(th)
    f1 = (ftest(th+h)-ftest(th-h))/(2*h)
    f2 = (ftest(th+h)-2*ftest(th)+ftest(th-h))/h**2
    f3 = (ftest(th+2*h)-2*ftest(th+h)+2*ftest(th-h)-ftest(th-2*h))/(2*h**3)
    f4 = (ftest(th+2*h)-4*ftest(th+h)+6*ftest(th)-4*ftest(th-h)+ftest(th-2*h))/h**4
    return f, f1, f2, f3, f4

def terms(th, f4):
    f, f1, f2, f3, _ = derivs(th)
    D, Dp = f2 + c0*f, f3 + c0*f1
    T, Tp, Tpp = -(1-p)*f1, -(1-p)*f2, -(1-p)*f3
    Dpp = f4 + c0*f2
    E = np.sqrt(0.75*D*D + 3*T*T)
    Ep = (1.5*D*Dp + 6*T*Tp)/(2*E)
    Epp = (1.5*(Dp**2 + D*Dpp) + 6*(Tp**2 + T*Tpp))/(2*E) - Ep**2/E
    G = E**(n-1)
    Gp = (n-1)*G*Ep/E
    Gpp = (n-1)*(Gp*Ep/E + G*Epp/E - G*Ep**2/E**2)
    return Gpp*D + 2*Gp*Dp + G*Dpp + a*G*D - 4*(Gp*T + G*Tp)   # should be 0

# f'''' from my rhs (which returns f4 = Dpp_rhs - c0*f2):
rhs = S.rhs_maker(p, n)
f, f1, f2, f3, _ = derivs(th0)
dy = rhs(th0, [f, f1, f2, f3])
f4_rhs = dy[3]

# f'''' required by the ODE, by root-finding on terms(th0, f4):
from scipy.optimize import brentq
g = lambda f4: terms(th0, f4)
lo, hi = f4_rhs - 50, f4_rhs + 50
f4_ode = brentq(g, lo, hi, xtol=1e-14)
print("f4_rhs  =", f4_rhs)
print("f4_ode  =", f4_ode)
print("rel diff =", abs(f4_rhs-f4_ode)/abs(f4_ode))
# FD 4th derivative as extra reference
_,_,_,_,f4_fd = derivs(th0)
print("f4_fd   =", f4_fd)

# also residual landscape near known HRR crack point n=4, alpha=pi, p=0.2
for q in [-0.5,-0.6,-0.7,-0.72,-0.75,-0.8,-0.9]:
    r = S.residual(0.2, q, 4.0, np.pi, rtol=1e-10, atol=1e-12)
    print(f"n=4 crack p=0.2 q={q:+.2f} -> f(a)={r[0]:+.3e} f'(a)={r[1]:+.3e}")
