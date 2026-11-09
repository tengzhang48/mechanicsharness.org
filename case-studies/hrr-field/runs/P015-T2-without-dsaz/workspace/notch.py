import numpy as np
from scipy.optimize import least_squares

# -----------------------------------------------------------------
# Leading plane-strain asymptotic field near a sharp V-notch tip in a
# power-law (Ramberg-Osgood, plastic-dominated) material.
#
#   sigma_ij = r^{-s} Sigma_ij(theta);   Airy  Phi = r^{2-s} f(theta)
#     Sigma_rr = f'' + (2-s) f
#     Sigma_tt = (2-s)(1-s) f
#     Sigma_rt = (s-1) f'
#   incompressible plastic:  eps_ij = (3/2) B sigma_e^{n-1} s_ij
#
# Reduced system (see notes):
#   P = C g^{n-1} D,   D = f'' + s(2-s) f,   R = (s-1) f'
#   g = sqrt( (3/4) D^2 + 3 R^2 )
#   U = (1-s n)^{-1} P
#   V' = -(2-s n)(1-s n)^{-1} P
#   U' - s n V - 4 C g^{n-1} R = 0
# C>0 is an overall factor; the eigenproblem is independent of it -> C=1.
# Unknowns: f_i, V_i on grid [0,alpha], eigenvalue s, gauge f_0=1.
# BCs: f'(0)=0, V(0)=0 (symmetry), Sigma_tt(alpha)=0, Sigma_rt(alpha)=0.
# -----------------------------------------------------------------

def solve_for(s_guess, n, alpha, N=400, f2_guess=0.0):
    th = np.linspace(0.0, alpha, N + 1)
    h = th[1] - th[0]

    # unknowns: f_0..f_N (N+1), V_0..V_N (N+1), s  -> but gauge f_0=1 fixed.
    # pack: x = [f_1..f_N (N), V_0..V_N (N+1), s]  ; f_0=1
    nf = N          # f_1..f_N
    nv = N + 1
    def unpack(x):
        f = np.empty(N + 1); f[0] = 1.0; f[1:] = x[:nf]
        V = x[nf:nf + nv]
        s = x[nf + nv]
        return f, V, s

    def fd1(a):
        out = np.empty_like(a)
        out[0] = (a[1] - a[0]) / h
        out[-1] = (a[-1] - a[-2]) / h
        out[1:-1] = (a[2:] - a[:-2]) / (2 * h)
        return out

    def fd2(a):  # second derivative, central interior, one-sided ends
        out = np.empty_like(a)
        out[1:-1] = (a[2:] - 2 * a[1:-1] + a[:-2]) / h**2
        out[0] = (a[2] - 2 * a[1] + a[0]) / h**2
        out[-1] = (a[-1] - 2 * a[-2] + a[-3]) / h**2
        return out

    def resid(x):
        f, V, s = unpack(x)
        D = fd2(f) + s * (2 - s) * f
        R = (s - 1) * fd1(f)
        g = np.sqrt(0.75 * D**2 + 3.0 * R**2)
        P = g**(n - 1) * D
        U = P / (1 - s * n)
        Vp = -(2 - s * n) / (1 - s * n) * P
        eps = 1e-9
        gs = np.where(g < eps, eps, g)  # avoid 0**0 issues
        RHS = 4.0 * gs**(n - 1) * R
        # eq1 (interior): V' - Vp = 0
        r1 = (fd1(V) - Vp)[1:-1] * h * 100
        # eq2 (interior): U' - s n V - RHS = 0
        r2 = (fd1(U) - s * n * V - RHS)[1:-1] * h * 100
        # BCs
        bc_fp0 = (fd1(f)[0]) * 1.0        # f'(0)=0
        bc_V0 = V[0] * 1.0                # V(0)=0
        # Sigma_tt(alpha)=0 ; Sigma_rt(alpha)=0
        Tt_a = (2 - s) * (1 - s) * f[-1]
        Rt_a = (s - 1) * fd1(f)[-1]
        return np.concatenate([r1, r2, [bc_fp0, bc_V0, Tt_a, Rt_a]])

    # initial guess
    f = np.cos(np.pi * th / (2 * alpha))  # smooth, f'(0)=0
    f = f / f[0]
    V = np.zeros(N + 1)
    x0 = np.concatenate([f[1:], V, [s_guess]])
    sol = least_squares(resid, x0, xtol=1e-13, ftol=1e-13, gtol=1e-13,
                        max_nfev=4000)
    f, V, s = unpack(sol.x)
    return s, f, V, th, sol

if __name__ == "__main__":
    n = 4.0
    alpha = 2 * np.pi / 3
    for sg in [0.1, 0.2, 0.3, 0.4, 0.5]:
        try:
            s, f, V, th, sol = solve_for(sg, n, alpha, N=300)
            print(f"s_guess={sg:5.2f} -> s={s:.6f}  cost={sol.cost:.3e}")
        except Exception as e:
            print("fail", sg, e)