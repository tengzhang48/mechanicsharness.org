"""Checks for the derivation panel of the poroelastic case page (6 Nov 2026).
Runs Codex's wilson_factor.py (unchanged, path given as argument) and checks, for nu = 0.2, nu_u = 0.4 (gamma = eta_t = 1/3):
the symbol N(xi) of VERDICT.md Eq. (4): N(0) = diag(1+gamma, 1), N(inf) = I, det N = 1 + gamma (m-k)/(m+k), Hermitian;
n11, n21 by three-grid extrapolation; K(0+)/sqrt(pi L) = n11/(1+gamma) (Eq. 6); n11^2 + n21^2 = 1 + gamma (Eq. 7);
A = A_faces + A_tips (Eqs. 9-11); the weak-coupling slope d(K/K_out)/d gamma -> -3/4 at gamma -> 0.
From the case package (NumPy); its output of 6 November 2026 is reference_checks/tip_derivation_check.txt:
    python3 reference_checks/tip_derivation_check.py grading/review_B_K10-codex/output/scratch/wilson_factor.py"""
import importlib.util, math, sys
import numpy as np
spec = importlib.util.spec_from_file_location('wf', sys.argv[1]); wf = importlib.util.module_from_spec(spec); spec.loader.exec_module(wf)
g = 1/3; nu, nuu = 0.2, 0.4; d = nuu - nu
assert abs(g - d/(1-nuu)) < 1e-15 and abs((1+g) - (1-nu)/(1-nuu)) < 1e-15
# the symbol, built as in wilson_factor.py lines 21-26, against Eq. (4)'s stated properties
xi = np.array([1e-9, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e6, -0.5, -7.0]); k = abs(xi); m = np.hypot(1, k)
N11 = 1 + g - 2*g*k*k/(m*(m+k)); N12 = np.sqrt(2*g)*k/(m+k)/(np.sqrt(1-1j*xi)*np.sqrt(-1j*xi)); N21 = N12.conj()
det = N11*1 - N12*N21
print('N11 at xi->0: %.12f (1+gamma = %.12f); at xi=1e6: %.9f (-> 1)' % (N11[0], 1+g, N11[6]))
print('|N12| at xi->0: %.2e; at xi=1e6: %.2e (both -> 0)' % (abs(N12[0]), abs(N12[6])))
print('max |det N - (1 + gamma (m-k)/(m+k))| = %.2e' % np.max(abs(det - (1 + g*(m-k)/(m+k)))))
print('min eigenvalue of N on the samples = %.6f (> 0: positive definite)' % min(np.linalg.eigvalsh(np.array([[N11[i], N12[i]], [N21[i], 1]])).min() for i in range(len(xi))))
rows = [wf.wilson(g, n, 0.01) for n in (16384, 65536, 262144)]
w = np.array([1/3, -2, 8/3]); n11 = float(w @ [r['n11'] for r in rows]); n21 = float(w @ [r['n21'] for r in rows])
for r in rows: print('count %6d: n11 = %.10f  n21 = %.10f  factor residual %.1e' % (r['count'], r['n11'], r['n21'], r['max_factor_residual']))
print('extrapolated: n11 = %.10f  n21 = %.10f' % (n11, n21))
K = n11/(1+g); print('K(0+)/sqrt(pi L) = n11/(1+gamma) = %.9f   (key 0.807593)' % K)
print('energy identity n11^2 + n21^2 - (1+gamma) = %.2e' % (n11**2 + n21**2 - (1+g)))
print('hydraulic share n21^2/(1+gamma) = %.4f%%' % (100*n21**2/(1+g)))
Af = 16*d/((1+g)*math.sqrt(math.pi)); At = 4*d/((1+g)*math.sqrt(math.pi))*math.pi*math.sqrt(2/g)*n21
print('A_faces = 16 d/((1+gamma) sqrt(pi)) = %.9f   (paper Eq. 18: 1.354055)' % Af)
print('A_tips  = 4 d/((1+gamma) sqrt(pi)) * pi sqrt(2/gamma) n21 = %.9f   (key 1.086165)' % At)
print('A = %.9f   (key 2.440220)' % (Af + At))
print('paper simplified K = 1/(1+gamma) = %.6f; elastic-energy estimate (1+gamma)^-1/2 = %.6f' % (1/(1+g), (1+g)**-0.5))
for e in (0.01, 0.02):
    r = [wf.wilson(e, n, 0.01) for n in (16384, 65536, 262144)]; n = float(w @ [x['n11'] for x in r])
    print('weak coupling: gamma = %.2f  K/K_out = %.6f  slope (K/K_out - 1)/gamma = %.4f  (expected -> -0.75)' % (e, n/(1+e), (n/(1+e) - 1)/e))
