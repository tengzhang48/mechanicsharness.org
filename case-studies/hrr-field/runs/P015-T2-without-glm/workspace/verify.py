import numpy as np
from scipy.optimize import least_squares
import solve as S

n, alpha = 4.0, 2*np.pi/3
p0, q0 = 0.1561752603, -1.26307900

# ---------- 1) tolerance/convergence of the shooting root ----------
print("tolerance study (shooting):")
for rt, at in [(1e-9,1e-11),(1e-10,1e-12),(1e-11,1e-13),(1e-12,1e-14),(1e-13,1e-15)]:
    s = least_squares(lambda x: S.residual(x[0], x[1], n, alpha, rt, at), [p0, q0],
                      xtol=1e-15, ftol=1e-15, gtol=1e-15)
    p, q = s.x
    print(f"  rtol={rt:.0e}: p={p:.10f} q={q:.8f} R={((2-p)+q)/((2-p)*(1-p)):.8f}")

# ---------- 2) independent collocation (cosine series) ----------
K, J = 16, 60
w = np.arange(K+1)*np.pi/(2*alpha)          # f'(0)=0 built in
thc = np.linspace(0.03, 0.97, J)*alpha      # interior collocation points

def basis(th, order):
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
    if np.any(E2 < 1e-10): return np.full(J+3, 1e6)
    E = np.sqrt(E2)
    Ep = (1.5*D*Dp + 6*T*Tp)/(2*E)
    Epp = (1.5*(Dp**2 + D*Dpp) + 6*(Tp**2 + T*Tpp))/(2*E) - Ep**2/E
    G = E**m; Gp = m*G*Ep/E
    Gpp = m*(Gp*Ep/E + G*Epp/E - G*Ep**2/E2)
    r = Gpp*D + 2*Gp*Dp + G*Dpp + aa*G*D - c1*(Gp*T + G*Tp)
    bc = [float((basis(np.array([alpha]), 0).T @ a)[0]),
          float((basis(np.array([alpha]), 1).T @ a)[0]),
          float((basis(np.array([0.0]), 0).T @ a)[0]) - 1.0]
    return np.concatenate([r, bc])

# initial guess from shooting solution
sol0 = S.shoot(p0, q0, n, alpha, rtol=1e-12, atol=1e-14)
ths = np.linspace(0, alpha, 801)
fs = sol0.sol(ths)[0]
a0 = np.linalg.lstsq(basis(ths, 0).T, fs, rcond=None)[0]
x0 = np.concatenate([a0, [p0]])
s = least_squares(resid_colloc, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
a, pcol = s.x[:K+1], s.x[K+1]
qcol = float((basis(np.array([0.0]), 2).T @ a)[0])
Rcol = ((2-pcol)+qcol)/((2-pcol)*(1-pcol))
print(f"collocation: p={pcol:.10f} q={qcol:.8f} R={Rcol:.8f} cost={s.cost:.3e} maxres={np.max(np.abs(s.fun)):.2e}")
print(f"vs shooting: p diff={abs(pcol-p0):.2e}  R diff={abs(Rcol-((2-p0)+q0)/((2-p0)*(1-p0))):.2e}")

# angular profile comparison (max abs diff of f)
fs_col = basis(ths, 0).T@a
print(f"max |f_shoot - f_colloc| = {np.max(np.abs(fs - fs_col)):.2e}")

# ---------- 3) summary quantities at theta=0 ----------
p, q = p0, q0
R0 = ((2-p)+q)/((2-p)*(1-p))
print(f"\nFINAL: p={p0:.7f} +-0.0003   R={R0:.6f} +-0.002")
print(f"  sigma_tt(0)~={(2-p)*(1-p):.5f} sigma_rr(0)~={2-p+q:.5f} sigma_zz(0)~={0.5*((2-p)+(2-p)*(1-p)):.5f}")
print(f"  strain exponent np = {n*p:.5f} (epsilon ~ r^-np), displacement r^(1-np) = r^{1-n*p:.5f}")
print(f"  energy integrability: p*(n+1) = {p*(n+1):.4f} < 2 OK;  J-like contour r^(1-p(n+1)) = r^{1-p*(n+1):.4f}")
