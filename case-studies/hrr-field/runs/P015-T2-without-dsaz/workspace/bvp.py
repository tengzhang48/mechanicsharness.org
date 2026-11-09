import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp

# =====================================================================
# Two-point BVP (collocation) for the plane-strain power-law V-notch
# leading mode-I field.  Unknown eigenvalue s.  Gauge f(0)=1.
# =====================================================================
s_sy, n_sy = sp.symbols('s n')
f0, f1, f2 = sp.symbols('f0 f1 f2')

D = f2 + s_sy*(2 - s_sy)*f0
R = (s_sy - 1)*f1
g = sp.sqrt(sp.Rational(3,4)*D**2 + 3*R**2)
P = g**(n_sy - 1)*D
GR = g**(n_sy - 1)*R
PF0, PF1, PF2 = sp.diff(P, f0), sp.diff(P, f1), sp.diff(P, f2)
GRF0, GRF1, GRF2 = sp.diff(GR, f0), sp.diff(GR, f1), sp.diff(GR, f2)
PF0_0, PF0_1, PF0_2 = sp.diff(PF0, f0), sp.diff(PF0, f1), sp.diff(PF0, f2)
PF1_0, PF1_1, PF1_2 = sp.diff(PF1, f0), sp.diff(PF1, f1), sp.diff(PF1, f2)
PF2_0, PF2_1 = sp.diff(PF2, f0), sp.diff(PF2, f1)
Vdot = -(2 - s_sy*n_sy)/(1 - s_sy*n_sy)*P
args = (f0, f1, f2, s_sy, n_sy)
L = lambda e: sp.lambdify(args, e, 'numpy')
PF0f, PF1f, PF2f = L(PF0), L(PF1), L(PF2)
GRF0f, GRF1f, GRF2f = L(GRF0), L(GRF1), L(GRF2)
PF0_0f, PF0_1f, PF0_2f = L(PF0_0), L(PF0_1), L(PF0_2)
PF1_0f, PF1_1f, PF1_2f = L(PF1_0), L(PF1_1), L(PF1_2)
PF2_0f, PF2_1f = L(PF2_0), L(PF2_1)
Pf, Vdotf = L(P), L(Vdot)

def make_fun(n):
    def fun(x, y, p):
        s = p[0]
        f0v, f1v, f2v, f3v, Vv = y
        a = (f0v, f1v, f2v, s, n)
        A = ((PF0_0f(*a)*f1v + PF0_1f(*a)*f2v + PF0_2f(*a)*f3v)*f1v
             + (PF1_0f(*a)*f1v + PF1_1f(*a)*f2v + PF1_2f(*a)*f3v)*f2v
             + (PF2_0f(*a)*f1v + PF2_1f(*a)*f2v)*f3v)
        B = PF2f(*a)
        Pv = Pf(*a)
        GRp = GRF0f(*a)*f1v + GRF1f(*a)*f2v + GRF2f(*a)*f3v
        f4 = (4.0*GRp - s*n*(2.0 - s)*Pv - A)/B
        Vp = Vdotf(*a)
        return np.vstack((f1v, f2v, f3v, f4, Vp))
    return fun

def make_bc(n):
    def bc(ya, yb, p):
        s = p[0]
        f0a, f1a, f2a, f3a, Va = ya
        f0b, f1b, f2b, f3b, Vb = yb
        Tt = (2.0 - s)*(1.0 - s)*f0b
        Rt = (s - 1.0)*f1b
        # normalize on curvature f''(0)=1 to exclude trivial (f=const, s=0) mode
        return np.array([f2a - 1.0, f1a, f3a, Va, Tt, Rt])
    return bc

def solve(n, alpha, s_guess, N=200):
    x = np.linspace(0, 1, N)
    th = alpha*x
    y0 = np.zeros((5, N))
    # seed: f''(0)=1 -> a curved even function; use cos(k th) with small decay
    y0[0] = np.cos(x*np.pi/2)
    y0[1] = -np.sin(x*np.pi/2)*np.pi/2
    y0[2] = -np.cos(x*np.pi/2)*(np.pi/2)**2
    y0[3] = np.sin(x*np.pi/2)*(np.pi/2)**3
    y0[4] = 0.0
    p0 = np.array([s_guess])
    sol = solve_bvp(make_fun(n), make_bc(n), th, y0, p=p0, tol=1e-9,
                    max_nodes=2000, verbose=0)
    return sol

if __name__ == "__main__":
    import time
    n = 4.0
    print("=== crack limit alpha->pi : expect s=1/(n+1)=0.2 ===")
    alpha = np.pi*0.999
    for sg in [0.2, 0.3, 0.5, 0.7]:
        t0=time.time()
        try:
            sol = solve(n, alpha, sg)
            print(f"  s_guess={sg}: success={sol.success} s={sol.p[0]:.6f} "
                  f"maxres={sol.rms_residuals.max():.2e} t={time.time()-t0:.1f}")
        except Exception as e:
            import traceback; traceback.print_exc()