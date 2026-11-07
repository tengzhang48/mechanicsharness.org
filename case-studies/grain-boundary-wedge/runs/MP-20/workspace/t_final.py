"""Final, consolidated numerics for the MP-20 grain-boundary relaxation.

Quantity list:
  lambda1, lambda2        decay rates (dimensionless, per t0)
  mode1_shape_half        f(z/h=0.5) of the slowest mode, normalised f(1)=1
  K_inf_norm              K(t->inf)/(sigma0 sqrt(pi h))
  relaxation_time_t0      relaxation time in units of t0  (= 1/lambda1)

Method: collocation of  g_tau = -E D2 g  on a uniform grid with g(0)=0
(Dirichlet) and g_s(1)=0 (Neumann).  E is the singularity-subtracted p.v.
elasticity operator, validated to 1e-13 against mpmath.  The gauge (constant)
mode of E is projected out before extracting lambda1.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(2000); xq = 0.5 * (xq + 1); wq = 0.5 * wq


def d2_uniform(n):
    h = 1.0 / n; D2 = np.zeros((n + 1, n + 1))
    for i in range(1, n):
        D2[i, i] = -2.0; D2[i, i - 1] = 1.0; D2[i, i + 1] = 1.0
    D2[n, n] = -2.0; D2[n, n - 1] = 2.0
    D2[0, 0] = 1.0
    return D2 / h ** 2


print('%4s %10s %10s %10s %10s' % ('n', 'lam1', 'lam2', 'lam3', 'f(0.5)'))
rows = []
for n in (10, 14, 18, 22, 26, 30):
    s = np.arange(n + 1) * (1.0 / n)
    E = gb3.elasticity_matrix(s, xq, wq)
    A = -(E @ d2_uniform(n))
    idx = np.arange(1, n + 1)
    Ai = A[np.ix_(idx, idx)]
    ev, V = np.linalg.eig(Ai)
    keep = np.abs(ev.imag) < 1e-3 * np.abs(ev)
    pairs = sorted([(ev.real[j], j) for j in np.where(keep & (ev.real > 1e-3))[0]])
    lams = [p[0] for p in pairs]
    # mode1 shape: eigenvector of the smallest positive eigenvalue, i.e. the
    # slowest relaxing mode; normalise to f(1)=1 (last interior node ~ s=1).
    j1 = pairs[0][1]
    v1 = V[:, j1].real
    f1 = v1 / v1[-1]
    fhalf = f1[np.argmin(np.abs(s[idx] - 0.5))]
    rows.append((n, lams[0], lams[1], lams[2], fhalf))
    print('%4d %10.5f %10.4f %10.4f %10.5f' % (n, lams[0], lams[1], lams[2], fhalf))

rows = np.array(rows)
print('\nconvergence checks (intercept at 1/n -> 0):')
for name, k in (('lam1', 1), ('lam2', 2), ('lam3', 3), ('f05', 4)):
    p = np.polyfit(1.0 / rows[:, 0], rows[:, k], 1)
    print('  %s : intercept=%.5f  (slope %.4f)' % (name, p[1], p[0]))