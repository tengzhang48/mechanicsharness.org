import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from notch import integrate_rk4, make_rhs, quadratic_coeff

n=4.0; alpha=2*np.pi/3
s=None; pf=None

# --- primary: RK4 shoot ---
x=fsolve(lambda x: np.array(integrate_rk4(n,x[0],x[1],alpha,nsteps=4000)),
         [1.8438247397,-1.26307900], xtol=1e-14)
s, pf = float(x[0]), float(x[1])

# --- independent: solve_ivp shoot to refine BC ---
def F_sv(x):
    sol=solve_ivp(lambda t,y:[y[1],y[2],y[3],make_rhs(n,x[0])(y[0],y[1],y[2],y[3])],
                  [0.0,alpha],[1.0,0.0,x[1],0.0],rtol=1e-13,atol=1e-15,max_step=1e-3)
    return np.array([sol.y[0,-1], sol.y[1,-1]])
x2=fsolve(F_sv,[s,pf],xtol=1e-13)
s2,pf2=float(x2[0]),float(x2[1])

print('RK4     s=%.12f fpp0=%.12f' % (s,pf))
print('solve_ivp s=%.12f fpp0=%.12f' % (s2,pf2))
print('delta_s=%.2e  delta_fpp0=%.2e' % (abs(s-s2), abs(pf-pf2)))

# stresses at theta=0 (Airy convention), f(0)=1, f'(0)=0, f''(0)=pf
sth0 = s*(s-1.0)*1.0
srr0 = s*1.0 + pf
R = srr0/sth0
p = 2.0 - s
print()
print('p = 2 - s = %.12f' % p)
print('sigma_theta(theta=0) = %.10f' % sth0)
print('sigma_r (theta=0)   = %.10f' % srr0)
print('R = srr0/sth0       = %.12f' % R)
