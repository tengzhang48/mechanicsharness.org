"""Mode-III (antiplane) crack-tip field for the power-law shear law of TASK.md.

gamma = (gam0/tau0) (tau/tau0)^(n-1) tau   =>   tau = tau0 (|grad w|/gam0)^(1/n)
Equilibrium div(tau)=0 => p-Laplacian div(|grad w|^{1/n-1} grad w)=0.
Separable leading field w = A r^s f(theta), s=1/(n+1); f solves the angular ODE and
f'(pi)=0 (traction-free faces), f odd in theta.  J fixes A.
"""
import numpy as np
from scipy.integrate import solve_ivp

n = 3.5; tau0 = 80.0; gam0 = 0.002
s = 1.0/(n+1.0); p = s; q = 1.0/n - 1.0

def _rhs(t, y):
    f, fp = y
    g2 = s*s*f*f + fp*fp
    return [fp, (-q*s*s*f*fp*fp - (1.0-p)*s*f*g2)/(g2 + q*fp*fp)]

sol = solve_ivp(_rhs, [0, np.pi], [0.0, 1.0], method='DOP853',
                rtol=1e-13, atol=1e-16, max_step=0.001, dense_output=True)

def ext(th):
    th = np.asarray(th, float)
    a = np.abs(th); f = sol.sol(a)[0]; fp = sol.sol(a)[1]
    return np.sign(th)*f, fp

def gext(th):
    f, fp = ext(th); return np.hypot(s*f, fp)

def ev(thr):
    a = np.array([abs(thr)]); f = sol.sol(a)[0][0]; fp = sol.sol(a)[1][0]
    if thr < 0: f = -f
    return f, fp, np.hypot(s*f, fp)

fpi = sol.sol(np.array([np.pi]))[0][0]

# J = C * A^{(n+1)/n} * I,  C = tau0*gam0^{-1/n}
_th = np.linspace(-np.pi, np.pi, 8000001)
_f, _fp = ext(_th); _g = gext(_th); _c = np.cos(_th); _sn = np.sin(_th)
_integ = ((n/(n+1.0))*_g**(1.0+1.0/n)*_c
          - _g**(1.0/n-1.0)*s*_f*(s*_f*_c - _fp*_sn))
I = np.trapezoid(_integ, _th)
C = tau0*gam0**(-1.0/n)
A = (0.120/(C*I))**(n/(n+1.0))

def tau_polar(r, thr):
    f, fp, g = ev(thr)
    pref = C*A**(1.0/n)*r**((s-1.0)/n)
    return pref*g**(1.0/n-1.0)*s*f, pref*g**(1.0/n-1.0)*fp

if __name__ == "__main__":
    print("n=%g s=1/(n+1)=%.10f  f(pi)=%.10f  I=%.10f  C=%.6f  A=%.12e"%(n,s,fpi,I,C,A))
    print("check f'(pi)=%.3e"%sol.sol(np.array([np.pi]))[1][0])
    # table
    print(" theta      f          f'         g")
    for d in (0,30,45,60,90,120,150,180):
        f, fp, g = ev(np.deg2rad(d))
        print("%6d  %9.5f  %9.5f  %9.5f"%(d,f,fp,g))
    print("Delta_w(0.0200)=%.10e mm   Delta_w(0.0800)=%.10e mm"%(2*A*0.02**s*fpi, 2*A*0.08**s*fpi))
    tr, tt = tau_polar(0.02, np.deg2rad(60.0))
    print("tau_rz(0.02,60)=%.8f MPa  tau_tz(0.02,60)=%.8f MPa"%(tr,tt))
    for d in (0,90,180):
        f, fp, g = ev(np.deg2rad(float(d)))
        print("r(|gamma|=0.01) theta=%3d : %.8f mm"%(d,(A*g/0.01)**(1.0/(1.0-s))))
    # direct Cartesian J verification (independent of the analytic coefficient)
    def direct_J(R, k=4000001):
        th = np.linspace(-np.pi, np.pi, k); f, fp = ext(th); g = gext(th)
        pref = C*A**(1.0/n)*R**((s-1.0)/n)
        tr = pref*g**(1.0/n-1.0)*(s*f); tt = pref*g**(1.0/n-1.0)*fp
        c = np.cos(th); sn = np.sin(th)
        W = (n/(n+1.0))*(pref*g**(1.0/n))*(A*R**(s-1)*g)
        wx = A*R**(s-1)*(s*f*c - fp*sn)
        return np.trapezoid((W*c - tr*wx)*R, th)
    print("direct J:", [round(direct_J(R), 9) for R in (0.02, 0.05, 0.1, 0.3)])
