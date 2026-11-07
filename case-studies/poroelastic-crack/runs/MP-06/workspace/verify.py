#!/usr/bin/env python3
"""Independent numerical verification of the analytic results in solve.py.

Parts:
  A. Chebyshev dislocation (crack) solver validated against closed-form
     pressurized-crack results (Sneddon): w(x), V, K, sigma_yy ahead of tip.
  B. Loading 2 steady state: 2D-FFT eigenstress pipeline; trace formula
     sigma_yy^sm(x,0) = -(gamma/2) p(x,0) for harmonic p; then K2_inf, V2_inf.
  C. Loading 1 transient: implicit-FD diffusion of p on half-plane with drained
     crack faces, eigenstress trace -> crack opening; V(t), K(t), short-time fit.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

nu, nu_u, alpha = 0.2, 0.4, 0.89
G, L = 1.0, 1.0
kappa = 3 - 4 * nu
gamma = alpha * (1 - 2 * nu) / (1 - nu)


# ---------------------------------------------------------------- crack solver
def crack_solve(q_nodes, Nc=241, Nm=160):
    """Solve PV Int mu(xi)/(x-xi) dxi = pi*(k+1)/(2G) q(x)  (uniform opening
    traction q on faces of crack |x|<L, traction-sign convention: total face
    traction to be cancelled is -q).  q_nodes given at Gauss-Chebyshev nodes.
    Returns theta, x, a_n, and callables for w, plus V, K."""
    j = np.arange(Nc)
    theta = (j + 0.5) * np.pi / Nc
    x = L * np.cos(theta)
    n = np.arange(1, Nm + 1)
    M = np.sin(np.outer(theta, n))            # M[j,m] = sin(n_m theta_j)
    b = -L * (kappa + 1) * q_nodes * np.sin(theta) / (2 * G)
    a, *_ = np.linalg.lstsq(M, b, rcond=None)
    # w(theta) = -sum a_n sin(n theta)/n
    def w_of_theta(th):
        return -np.sum(a[None, :] * np.sin(np.outer(th, n)) / n[None, :], axis=1)
    # V = -(pi L/2) sum_{n odd} a_n / n
    odd = (n % 2 == 1)
    V = -np.pi * L / 2 * np.sum(a[odd] / n[odd])
    # K = (1/sqrt(pi L)) Int q(x) sqrt((L+x)/(L-x)) dx  (Gauss-Chebyshev quad)
    K = np.sum(q_nodes * np.sqrt((L + x) / (L - x)) * np.pi / Nc * np.sin(theta)) / np.sqrt(np.pi * L)
    return theta, x, a, w_of_theta, V, K, n


def yy_ahead_of_tip(a, n, xa):
    """sigma_yy(x,0) for x>L from dislocation density (for validation)."""
    # PV Int mu/(x-xi) with mu = F/sqrt(L^2-xi^2), F = sum a_n cos(n phi):
    # Int_0^pi cos(n phi)/(x/L - cos phi) dphi = pi (r^-n - r^n)/( ... ) classic:
    # for X=x/L>1: PV Int_0^pi cos(n phi)/(X-cos phi) dphi = pi*( (X-sqrt(X^2-1))^n )*(1)*? sign
    # Standard: PV Int cos(n phi)/(X-cos phi) dphi = -pi*(X - sqrt(X^2-1))^n/sqrt(X^2-1)*(X... )
    # use numerical principal-value quadrature instead (robust):
    Nq = 40000
    phi = (np.arange(Nq) + 0.5) * np.pi / Nq
    xi = L * np.cos(phi)
    F = np.cos(np.outer(phi, n)) @ a
    # PV Int mu(xi)/(x-xi) dxi = Int_0^pi F(phi)/(x - L cos(phi)) dphi
    out = []
    for X in np.atleast_1d(xa):
        integrand = F / (X - L * np.cos(phi))
        out.append(np.sum(integrand) * np.pi / Nq)
    Pv = np.array(out)
    return -2 * G / (np.pi * (kappa + 1)) * Pv


print("=" * 78)
print("PART A: crack (dislocation) solver validation, uniform face traction q=1")
print("=" * 78)
Nc = 241
th, xq, a, w_of_theta, V, K, nn = crack_solve(np.ones(Nc))
thx = np.linspace(0.05, np.pi - 0.05, 13)
xx = L * np.cos(thx)
w_num = w_of_theta(thx)
w_ex = 2 * (1 - nu) * np.sqrt(L * L - xx * xx) / G
print("max |w_num - w_exact| / max|w_exact| = %.2e"
      % (np.max(np.abs(w_num - w_ex)) / np.max(np.abs(w_ex))))
print("V  = %.6f   exact pi*(1-nu) = %.6f   (err %.2e)" % (V, np.pi * (1 - nu), abs(V - np.pi * (1 - nu))))
print("K  = %.6f   K/sqrt(pi L) = %.6f  (exact 1.000000, err %.2e)"
      % (K, K / np.sqrt(np.pi * L), abs(K / np.sqrt(np.pi * L) - 1.0)))
xa = np.array([1.05, 1.2, 1.5, 2.0, 3.0])
syy = yy_ahead_of_tip(a, nn, xa)
syy_ex = xa / np.sqrt(xa * xa - L * L) - 1.0     # Westergaard pressurized crack
print("sigma_yy(x,0), x>L  numeric vs Westergaard closed form:")
for X, s1, s2 in zip(xa, syy, syy_ex):
    print("   x=%.2f  num=%9.5f  exact=%9.5f  reldiff=%.2e" % (X, s1, s2, abs(s1 - s2) / abs(s2)))

# ---------------------------------------------------------------- part B
print()
print("=" * 78)
print("PART B: loading 2 steady state (p_ss=1 on faces, harmonic)")
print("=" * 78)
A_box, Nf = 4.0, 2048
xs = (np.arange(Nf) - Nf // 2) * (2 * A_box / Nf)
ys = xs.copy()
Xg, Yg = np.meshgrid(xs, ys, indexing="xy")
P = (np.arctan2(Xg + L, np.abs(Yg)) + np.arctan2(L - Xg, np.abs(Yg))) / np.pi
kx = 2 * np.pi * np.fft.fftfreq(Nf, d=2 * A_box / Nf)
ky = kx.copy()
KX, KY = np.meshgrid(kx, ky, indexing="xy")
K2 = KX**2 + KY**2
ker = np.zeros_like(K2)
ker[K2 > 0] = KX[K2 > 0] ** 2 / K2[K2 > 0]   # kx^2/(kx^2+ky^2)
Ph = np.fft.fft2(P)
Syy = np.fft.ifft2(-gamma * ker * Ph).real
row = Nf // 2        # y = 0
sel = np.abs(xs) < L - 0.02
num_trace = Syy[row, :]
an_trace = -gamma / 2 * P[row, :]
err = np.max(np.abs(num_trace[sel] - an_trace[sel])) / gamma
print("FFT trace sigma_yy^sm(x,0) vs -(gamma/2) p(x,0) inside faces (away from tips):")
print("   max rel err = %.3f   (gamma=%.4f, predicted value -gamma/2=%.4f)"
      % (err, gamma, -gamma / 2))
print("   sample: p_ss(x=0.5)=%.4f, FFT trace=%.4f, formula=%.4f"
      % (P[row, np.argmin(np.abs(xs - 0.5))], num_trace[np.argmin(np.abs(xs - 0.5))],
         an_trace[np.argmin(np.abs(xs - 0.5))]))
# crack solution with q = -sigma_yy^sm = gamma/2 (uniform): from PART A linearity
thq, xq2, _, _, V2, K2v, _ = crack_solve(np.full(Nc, -gamma / 2))
print("K2_inf  crack-solver (q=-gamma/2) = %.6f /sqrt(pi L) -> %.6f   analytic -gamma/2 = %.6f"
      % (K2v, K2v / np.sqrt(np.pi * L), -gamma / 2))
print("V2_inf  crack-solver (q=-gamma/2) = %.6f   analytic -pi*alpha*(1-2nu)/2 = %.6f"
      % (V2, -np.pi * alpha * (1 - 2 * nu) / 2))
# same but with q taken from the FFT trace directly
qF = num_trace   # solver input is sigma_yy^sm itself (negative on faces)
q_nodes = np.interp(xq, xs, qF)
_, _, _, _, V2F, K2F, _ = crack_solve(q_nodes)
print("K2_inf  with FFT-derived q      = %.6f /sqrt(pi L) -> %.6f" % (K2F, K2F / np.sqrt(np.pi * L)))
print("V2_inf  with FFT-derived q      = %.6f" % V2F)
