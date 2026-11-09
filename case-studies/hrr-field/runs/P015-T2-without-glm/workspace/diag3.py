import numpy as np
import solve as S

p, q, n, alpha = 0.1561752603, -1.26307900, 4.0, 2*np.pi/3
sol = S.shoot(p, q, n, alpha, rtol=1e-13, atol=1e-15)

# A) pointwise: does dense-output trajectory satisfy the same f'''' the integrator used?
th = np.linspace(0.05*alpha, 0.95*alpha, 4001)
Y = sol.sol(th)
f4_fd = np.gradient(Y[3], th)  # d f'''/d theta
rhs4 = np.array([S.rhs_maker(p, n)(t, Y[:, i])[3] for i, t in enumerate(th)])
i = np.argmax(np.abs(f4_fd - rhs4))
j = np.argsort(np.abs(f4_fd - rhs4))[len(th)//2]
for k in [j, i]:
    print(f"th={th[k]:.4f} f4_fd={f4_fd[k]:+.6f} rhs_f4={rhs4[k]:+.6f} rel={abs(f4_fd[k]-rhs4[k])/max(1e-9,abs(rhs4[k])):.3e}")

# B) verify rhs f'''' against the NEW ODE on a smooth test function
n_, p_ = 4.0, 0.1561752603
c0 = p_*(2-p_); a = n_*p_*(2-n_*p_); c1 = 4*(1-n_*p_)
def ftest(th): return 1.0 + 0.3*np.cos(1.7*th) - 0.2*np.cos(0.9*th)**2
h = 1e-3
def d(th, k):  # k-th derivative via central FD
    st = {1:[-1,1],2:[1,-2,1],3:[-1,2,-2,1],4:[1,-4,6,-4,1]}[k]
    c = {1:[-0.5,0.5],2:[1,-2,1],3:[-0.5,1,-1,0.5],4:[1,-4,6,-4,1]}[k]
    return sum(cc*ftest(th+ii*h) for cc, ii in zip(c, {1:[-1,1],2:[-1,0,1],3:[-2,-1,1,2],4:[-2,-1,0,1,2]}[k]))/h**k
th0 = 0.7
f, f1, f2, f3 = ftest(th0), d(th0,1), d(th0,2), d(th0,3)
D, Dp = f2 + c0*f, f3 + c0*f1
T, Tp, Tpp = -(1-p_)*f1, -(1-p_)*f2, -(1-p_)*f3
E = np.sqrt(0.75*D*D + 3*T*T)
Ep = (1.5*D*Dp + 6*T*Tp)/(2*E)
G = E**(n_-1)
Gp = (n_-1)*G*Ep/E
# ODE as function of f4:  G''D + 2G'D' + GD'' + aGD - c1(G'T+GT') = 0
def terms(f4):
    Dpp = f4 + c0*f2
    Epp = (1.5*(Dp**2 + D*Dpp) + 6*(Tp**2 + T*Tpp))/(2*E) - Ep**2/E
    Gpp = (n_-1)*(Gp*Ep/E + G*Epp/E - G*Ep**2/E**2)
    return Gpp*D + 2*Gp*Dp + G*Dpp + a*G*D - c1*(Gp*T + G*Tp)
from scipy.optimize import brentq
f4_rhs = S.rhs_maker(p_, n_)(th0, [f, f1, f2, f3])[3]
f4_ode = brentq(terms, f4_rhs-100, f4_rhs+100, xtol=1e-13)
print(f"B) f4_rhs={f4_rhs:.10f} f4_ode={f4_ode:.10f} rel={abs(f4_rhs-f4_ode)/abs(f4_ode):.2e}")

# C) same but on the trajectory point used in A
Y0 = sol.sol(np.array([th0]))[:, 0]
f, f1, f2, f3 = Y0
D, Dp = f2 + c0*f, f3 + c0*f1
T, Tp, Tpp = -(1-p_)*f1, -(1-p_)*f2, -(1-p_)*f3
E = np.sqrt(0.75*D*D + 3*T*T)
Ep = (1.5*D*Dp + 6*T*Tp)/(2*E)
G = E**(n_-1); Gp = (n_-1)*G*Ep/E
f4_rhs = S.rhs_maker(p_, n_)(th0, Y0)[3]
f4_ode = brentq(terms, f4_rhs-100, f4_rhs+100, xtol=1e-13)
print(f"C) f4_rhs={f4_rhs:.10f} f4_ode={f4_ode:.10f} rel={abs(f4_rhs-f4_ode)/abs(f4_ode):.2e}")

# D) FD residual of full composite along trajectory, interior only, relative to local scale
c0 = p*(2-p); a = n*p*(2-n*p); c1 = 4*(1-n*p)
D = Y[2] + c0*Y[0]; T = -(1-p)*Y[1]
E = np.sqrt(0.75*D*D + 3*T*T); G = E**(n-1)
GD = G*D; GT = G*T
GDpp = np.gradient(np.gradient(GD, th), th)
GTp = np.gradient(GT, th)
L = GDpp + a*GD - c1*GTp
loc = np.maximum(np.abs(a*GD), np.maximum(np.abs(c1*GTp), np.abs(GDpp)))
loc[loc < 1e-2*np.max(loc)] = np.max(loc)
print(f"D) interior ODE rel resid: max={np.max(np.abs(L)/loc):.3e} median={np.median(np.abs(L)/loc):.3e}")
ii = np.argmax(np.abs(L)/loc)
print(f"   worst at th={th[ii]:.4f} (alpha={alpha:.4f}); L={L[ii]:.4e} terms: GDpp={GDpp[ii]:.4e} aGD={a*GD[ii]:.4e} c1GTp={c1*GTp[ii]:.4e}")
