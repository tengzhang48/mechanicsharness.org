import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# Material / loading
tau0 = 80.0
gamma0 = 0.002
n = 3.5
J = 0.120

p = 1.0/(n+1.0)   # stress exponent
s = 1.0/(n+1.0)   # displacement exponent (antisymmetric/antiplane HRR)
q = 1.0/n         # tau ~ gamma^(1/n)

# Angular ODE.  write tau = C A^(1/n) r^(-p) T(theta)
# T_r = g^(1/n-1) * (s f),  T_th = g^(1/n-1) * f',  g = sqrt((s f)^2 + f'^2)
# equilibrium: (1-p) T_r + T_th' = 0
def rhs(theta, y):
    f, fp = y
    g2 = (s*f)**2 + fp**2
    g = np.sqrt(max(g2, 1e-300))
    gm = g**(1.0/n - 1.0)
    Tr = gm * (s*f)
    dTth_df = gm * (1.0/n - 1.0) * (s**2 * f) / g2 * fp   # d/dtheta via df; not needed if we use dTth/dtheta
    # T_th' = d/dtheta[ g^(1/n-1) f' ] ; compute directly
    # dg/dtheta = (s^2 f fp + fp f'')/g ; but f'' from ODE. Instead differentiate expression symbolically:
    # T_th' = (1/n-1) g^(1/n-2) g' fp + g^(1/n-1) f''
    # and g g' = (s^2 f fp + fp f'')  -> contains f''. Solve for f'':
    # Let k = 1/n - 1. T_th' = g^k f'' + k g^(k-1) g' fp
    # g' = (s^2 f fp + fp f'')/g
    # T_th' = g^k f'' + k g^(k-2) fp (s^2 f fp + fp f'')
    #       = f'' [g^k + k g^(k-2) fp^2] + k g^(k-2) s^2 f fp^2
    # equilibrium: (1-p) g^k s f + f''[g^k + k g^(k-2) fp^2] + k g^(k-2) s^2 f fp^2 = 0
    k = 1.0/n - 1.0
    A22 = g**k + k * g**(k-2.0) * fp**2
    b = (1.0-p) * g**k * s * f + k * g**(k-2.0) * s**2 * f * fp**2
    fpp = -b / A22
    return [fp, fpp]

deg = np.pi/180.0
sol = solve_ivp(rhs, [0.0, np.pi], [0.0, 1.0], rtol=1e-12, atol=1e-14, dense_output=True, method='DOP853')
f_pi = sol.y[0, -1]
fp_pi = sol.y[1, -1]
print("f(pi)=", f_pi, " f'(pi)=", fp_pi, " (should be ~0)")
print("f(pi/2)=", sol.sol(np.pi/2)[0], " f'(pi/2)=", sol.sol(np.pi/2)[1])

# g at key angles
def gfun(th):
    f, fp = sol.sol(th)
    return np.sqrt((s*f)**2 + fp**2)
for th in [0.0, np.pi/2, np.pi]:
    print("g(", th/deg, "deg)=", gfun(th))
