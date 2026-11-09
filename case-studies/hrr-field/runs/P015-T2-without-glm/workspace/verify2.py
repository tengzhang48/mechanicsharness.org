import numpy as np
from scipy.optimize import least_squares
import solve as S

n, alpha = 4.0, 2*np.pi/3
p0, q0 = 0.1561752603, -1.26307900
R0 = ((2-p0)+q0)/((2-p0)*(1-p0))

# ---------- A) collocation, fixed scaling ----------
K, J = 12, 48
w = np.arange(K+1)*np.pi/(2*alpha)
thc = np.linspace(0.03, 0.97, J)*alpha

def basis(th, order):
    th = np.atleast_1d(th)
    if order == 0: return np.cos(np.outer(w, th))
    if order == 1: return -np.sin(np.outer(w, th))*(w[:, None])
    if order == 2: return -np.cos(np.outer(w, th))*(w[:, None]**2)
    if order == 3: return np.sin(np.outer(w, th))*(w[:, None]**3)
    if order == 4: return np.cos(np.outer(w, th))*(w[:, None]**4)

def resid_colloc(x):
    a = x[:K+1]; p = x[K+1]
    c0 = p*(2-p); aa = n*p*(2-n*p); c1 = 4*(1-n*p); m = n-1
    F0 = basis(thc,0).T@a; F1 = basis(thc,1).T@a; F2 = basis(thc,2).T@a
    F3 = basis(thc,3).T@a; F4 = basis(thc,4).T@a
    D, Dp, Dpp = F2+c0*F0, F3+c0*F1, F4+c0*F2
    T, Tp, Tpp = -(1-p)*F1, -(1-p)*F2, -(1-p)*F3
    E2 = 0.75*D*D + 3*T*T
    if np.any(E2 < 1e-9): return np.full(J+3, 1e6)
    E = np.sqrt(E2)
    Ep = (1.5*D*Dp + 6*T*Tp)/(2*E)
    Epp = (1.5*(Dp**2 + D*Dpp) + 6*(Tp**2 + T*Tpp))/(2*E) - Ep**2/E
    G = E**m; Gp = m*G*Ep/E
    Gpp = m*(Gp*Ep/E + G*Epp/E - G*Ep**2/E2)
    r = Gpp*D + 2*Gp*Dp + G*Dpp + aa*G*D - c1*(Gp*T + G*Tp)
    bc = [float((basis(alpha,0).T@a)[0]), float((basis(alpha,1).T@a)[0]),
          float((basis(0.0,0).T@a)[0]) - 1.0]
    return np.concatenate([r, bc])

sol0 = S.shoot(p0, q0, n, alpha, rtol=1e-12, atol=1e-14)
ths = np.linspace(0, alpha, 2001)
fs = sol0.sol(ths)[0]
a0 = np.linalg.lstsq(basis(ths, 0).T, fs, rcond=None)[0]
x0 = np.concatenate([a0, [p0]])
r0v = resid_colloc(x0)
print("colloc initial: max|res_ode|=%.3e max|res_bc|=%.3e  fit_err=%.2e"
      % (np.max(np.abs(r0v[:J])), np.max(np.abs(r0v[J:])), np.max(np.abs(basis(ths,0).T@a0 - fs))))
s = least_squares(resid_colloc, x0, x_scale='jac', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=40000)
a, pcol = s.x[:K+1], s.x[K+1]
qcol = float((basis(0.0, 2).T @ a)[0])
Rcol = ((2-pcol)+qcol)/((2-pcol)*(1-pcol))
print(f"collocation: p={pcol:.10f} R={Rcol:.8f} cost={s.cost:.3e} maxres={np.max(np.abs(s.fun)):.2e}")
print(f"vs shooting : p diff={abs(pcol-p0):.2e}  R diff={abs(Rcol-R0):.2e}")

# ---------- B) flank-side shooting (independent IVP direction) ----------
# y(alpha) = [0,0,s,t];  residuals: f'(0)=0, f'''(0)=0, f(0)=1  -> unknowns s,t,p
from scipy.integrate import solve_ivp
def rhs_neg(p):
    base = S.rhs_maker(p, n)
    def rhs(th, y):
        return [-v for v in base(-th, y)]
    return rhs
def resid_flank(x):
    s, t, p = x
    try:
        sol = solve_ivp(rhs_neg(p), (0.0, alpha), [0.0, 0.0, s, t],
                        method='DOP853', rtol=1e-12, atol=1e-14)
    except Exception:
        return np.array([1e6]*3)
    if not sol.success: return np.array([1e6]*3)
    f0, f1_0, _, f3_0 = sol.y[:, -1]   # values at theta=0
    return np.array([f1_0, f3_0, f0 - 1.0])

# seed: f''(a), f'''(a) from dense output slightly inside the flank
tt = np.linspace(alpha-1e-4, alpha, 3)
Yt = sol0.sol(tt)
sf, tf = float(Yt[2, 0]), float(np.gradient(Yt[3], tt)[0])
s2 = least_squares(resid_flank, [sf, tf, p0], x_scale='jac', xtol=1e-15, ftol=1e-15, gtol=1e-15)
s_, t_, pf = s2.x
# R from flank solution: f(0)=1 => same normalization; sigma_rr(0) = (2-p)f(0)+f''(0), need f''(0):
solb = solve_ivp(rhs_neg(pf), (0.0, alpha), [0.0, 0.0, s_, t_], method='DOP853', rtol=1e-12, atol=1e-14)
qf = float(solb.y[2, -1])
Rf = ((2-pf)+qf)/((2-pf)*(1-pf))
print(f"flank shoot : p={pf:.10f} R={Rf:.8f} res={s2.fun} cost={s2.cost:.2e}")
print(f"vs shooting : p diff={abs(pf-p0):.2e}  R diff={abs(Rf-R0):.2e}")
