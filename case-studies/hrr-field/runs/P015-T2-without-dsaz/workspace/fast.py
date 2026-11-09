import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

s_sy, n_sy = sp.symbols('s n')
f0, f1, f2 = sp.symbols('f0 f1 f2')

D = f2 + s_sy*(2 - s_sy)*f0
R = (s_sy - 1)*f1
g = sp.sqrt(sp.Rational(3,4)*D**2 + 3*R**2)
P = g**(n_sy - 1)*D
GR = g**(n_sy - 1)*R

# partials of P and GR wrt (f0,f1,f2)
PF0 = sp.diff(P, f0); PF1 = sp.diff(P, f1); PF2 = sp.diff(P, f2)
GRF0 = sp.diff(GR, f0); GRF1 = sp.diff(GR, f1); GRF2 = sp.diff(GR, f2)
# Second-level partials needed for P'': P'' = d/dth (PF0 f1 + PF1 f2 + PF2 f3)
# Only PF2 (=dP/df2) multiplies f4 (since df3/dth = f4 and PF2 depends on f0,f1,f2).
# d/dth PFj = (dPFj/df0) f1 + (dPFj/df1) f2 + (dPFj/df2) f3
PF0_0 = sp.diff(PF0, f0); PF0_1 = sp.diff(PF0, f1); PF0_2 = sp.diff(PF0, f2)
PF1_0 = sp.diff(PF1, f0); PF1_1 = sp.diff(PF1, f1); PF1_2 = sp.diff(PF1, f2)
PF2_0 = sp.diff(PF2, f0); PF2_1 = sp.diff(PF2, f1)  # PF2_2 will pair with f4

# V-dot: -(2 - s n)/(1 - s n) * P
Vdot = -(2 - s_sy*n_sy)/(1 - s_sy*n_sy)*P

args = (f0, f1, f2, s_sy, n_sy)
def L(expr):
    return sp.lambdify(args, expr, 'numpy')

PF0f, PF1f, PF2f = L(PF0), L(PF1), L(PF2)
GRF0f, GRF1f, GRF2f = L(GRF0), L(GRF1), L(GRF2)
PF0_0f, PF0_1f, PF0_2f = L(PF0_0), L(PF0_1), L(PF0_2)
PF1_0f, PF1_1f, PF1_2f = L(PF1_0), L(PF1_1), L(PF1_2)
PF2_0f, PF2_1f = L(PF2_0), L(PF2_1)
Pf = L(P)
Vdotf = L(Vdot)

def make_rhs(n, s):
    def rhs(th, y):
        f0v, f1v, f2v, f3v, Vv = y
        a = (f0v, f1v, f2v, s, n)
        # P'' = A + B*f4
        # P' = PF0*f1 + PF1*f2 + PF2*f3 ; d/dth:
        A = ((PF0_0f(*a)*f1v + PF0_1f(*a)*f2v + PF0_2f(*a)*f3v)*f1v
             + (PF1_0f(*a)*f1v + PF1_1f(*a)*f2v + PF1_2f(*a)*f3v)*f2v
             + (PF2_0f(*a)*f1v + PF2_1f(*a)*f2v)*f3v)
        B = PF2f(*a)
        Pv = Pf(*a)
        GRp = GRF0f(*a)*f1v + GRF1f(*a)*f2v + GRF2f(*a)*f3v
        # ODE: P'' + s n (2-s) P - 4 GR' = 0
        f4 = (4.0*GRp - s*n*(2.0 - s)*Pv - A)/B
        Vp = Vdotf(*a)
        return [f1v, f2v, f3v, f4, Vp]
    return rhs

def integrate(s, n, alpha, f2_0, f0_0=1.0):
    # mode I symmetric: f even -> f'(0)=0 (=f1) and f'''(0)=0 (=f3); V odd -> V(0)=0
    y0 = [f0_0, 0.0, f2_0, 0.0, 0.0]
    return solve_ivp(make_rhs(n, s), [0.0, alpha], y0, rtol=1e-9,
                     atol=1e-11, method='RK45', dense_output=True, max_step=0.05)

def bcs(x, n, alpha):
    s, f2_0 = x
    sol = integrate(s, n, alpha, f2_0)
    f0v, f1v, f2v, f3v, Vv = sol.y[:, -1]
    Tt = (2.0 - s)*(1.0 - s)*f0v     # Sigma_tt(alpha)=0
    Rt = (s - 1.0)*f1v               # Sigma_rt(alpha)=0
    return np.array([Tt, Rt])

def solve(n, alpha, s0, f2_0=1.0):
    from scipy.optimize import root
    r = root(bcs, [s0, f2_0], args=(n, alpha), method='hybr',
             options={'xtol': 1e-13, 'maxfev': 800})
    return r.x, r.fun, r.success

if __name__ == "__main__":
    import time
    t0=time.time()
    alpha = np.pi*0.999
    print("crack limit n=4 target 0.2")
    roots=set()
    for s0 in np.linspace(0.02,0.95,30):
        try:
            x,r,ok=solve(4.0,alpha,s0,f2_0=1.0)
            if np.linalg.norm(r)<1e-9 and x[0]>1e-3:
                roots.add(round(x[0],4))
        except Exception: pass
    print("roots:",sorted(roots),"time",round(time.time()-t0,1))