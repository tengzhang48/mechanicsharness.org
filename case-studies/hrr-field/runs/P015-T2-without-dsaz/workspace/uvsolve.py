import numpy as np
from scipy.optimize import least_squares

# Incompressible power-law plane-strain V-notch, DISPLACEMENT formulation.
# Unknowns on theta in [0, alpha]: U, V, P at M nodes, plus eigenvalue s.
# Equations (shapes, B=1):
#   dev: e2 = (2/3)(d_rr^2 + d_tt^2 + 2 ert^2),  K = e2^{(1-n)/(2n)}
#   s_ij = (2/3) K * eps_ij ;  S_ij = s_ij + P delta
#   eqr: -s n S_rr + (S_rr - S_tt) + S_rt' = 0
#   eqt: (-s n + 2) S_rt + S_tt' = 0
#   inc: (2 - s n) U + V' = 0
# with
#   d_rr = (1 - s n) U
#   d_tt = U + V'
#   ert  = 1/2 (U' - s n V)
# Hmm: use inc to eliminate. Actually keep d_tt and inc consistent.

def make(n, alpha, M):
    th = np.linspace(0, alpha, M)
    h = th[1]-th[0]
    def D1(y):  # first derivative, 2nd order interior
        d = np.zeros_like(y)
        d[1:-1] = (y[2:]-y[:-2])/(2*h)
        d[0] = (y[1]-y[0])/h
        d[-1] = (y[-1]-y[-2])/h
        return d
    def D2(y):
        d = np.zeros_like(y)
        d[1:-1] = (y[2:]-2*y[1:-1]+y[:-2])/h**2
        d[0] = (y[2]-2*y[1]+y[0])/h**2
        d[-1] = (y[-2]-2*y[-1]+y[-1])/h**2
        return d
    return th, D1, D2

def residual(z, n, alpha, M, sfix=None):
    th, D1, D2 = make(n, alpha, M)
    U = z[0:M]; V = z[M:2*M]; P = z[2*M:3*M]
    s = z[3*M] if sfix is None else sfix
    Up = D1(U); Vp = D1(V); Upp = D2(U); Vpp = D2(V); Pp = D1(P)
    d_rr = (1 - s*n)*U
    d_tt = U + Vp
    ert  = 0.5*(Up - s*n*V)
    e2 = (2/3.0)*(d_rr**2 + d_tt**2 + 2*ert**2) + 1e-30
    K = e2**((1-n)/(2.0*n))
    s_rr = (2/3.0)*K*d_rr
    s_tt = (2/3.0)*K*d_tt
    s_rt = (2/3.0)*K*ert
    S_rr = s_rr + P; S_tt = s_tt + P; S_rt = s_rt
    S_rt_p = D1(S_rt); S_tt_p = D1(S_tt)
    eqr = -s*n*S_rr + (S_rr - S_tt) + S_rt_p
    eqt = (-s*n + 2)*S_rt + S_tt_p
    inc = (2 - s*n)*U + Vp
    # ---- boundary / symmetry conditions -------------------------------
    # theta=0 (mode I symmetry):  V(0)=0 ; U'(0)=0 ; (also s_rt(0)=0 auto)
    # theta=alpha (traction free): S_tt(alpha)=0 ; S_rt(alpha)=0
    rows = np.concatenate([eqr, eqt, inc])
    # replace a few rows with BCs (keep problem square-ish; use least squares)
    # gauge: fix amplitude -> U[0] = 1
    bcs = np.array([
        V[0],                         # symmetry flux
        Up[0],                        # symmetry
        S_tt[-1],                     # traction
        S_rt[-1],                     # traction
        U[0] - 1.0,                   # amplitude gauge
    ])
    # weight BCs strongly, and drop the constraint rows at the boundary ends
    full = np.concatenate([rows, bcs*10.0])
    return full

def solve(n, alpha, s0, M=201, sfix=None):
    z0 = np.zeros(3*M+1)
    th = np.linspace(0, alpha, M)
    # seed from elastic-like: U ~ cos, V ~ sin
    z0[0:M] = np.cos(1.5*th)
    z0[M:2*M] = -np.sin(th)
    z0[2*M:3*M] = 0.0
    z0[3*M] = s0
    kw = dict(method='lm', max_nfev=20000, xtol=1e-13, ftol=1e-13)
    if sfix is not None:
        fun = lambda z: residual(z, n, alpha, M, sfix=sfix)
    else:
        fun = lambda z: residual(z, n, alpha, M)
    sol = least_squares(fun, z0, **kw)
    return sol

def solve2(n, alpha, s0, M=81):
    th = np.linspace(0, alpha, M)
    z0 = np.zeros(3*M+1)
    z0[0:M] = np.cos(1.5*th)
    z0[M:2*M] = -np.sin(th)
    z0[3*M] = s0
    sol = least_squares(residual, z0, args=(n, alpha, M),
                        method='lm', max_nfev=4000, xtol=1e-12, ftol=1e-12)
    return sol

if __name__ == "__main__":
    import sys
    n = float(sys.argv[1]); a = float(sys.argv[2]); s0 = float(sys.argv[3]); M = int(sys.argv[4])
    sol = solve2(n, a, s0, M=M)
    s = sol.x[3*M]
    print(f"n={n} alpha={a:.4f} M={M} -> s={s:.6f} cost={sol.cost:.3e} res={np.linalg.norm(sol.fun):.3e}")