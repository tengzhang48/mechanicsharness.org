import numpy as np
from scipy.integrate import quad
from sblib import M_moments, G_moments, R_all, c_coeff, build_s_matrices, f_from_coeffs
from gblib import S_kernel, T_moment

print("=== S1: R_m vs direct PV quad of full kernel ===")
def R_quad(s, m):
    x = 1.0-s*s
    # PV at r=s: subtract r^m/(r^2-s^2), add back M_m analytic
    def integrand(r):
        xi = 1.0-r*r
        return r**m*(S_kernel(x, xi) - 1.0/(x-xi))
    val, _ = quad(integrand, 0, 1, epsabs=1e-12, epsrel=1e-12, limit=300)
    M = M_moments(s, m)[m]
    return val + M

for s in [0.05, 0.3, 0.65, 0.9, 0.99, 0.99999]:
    for m in [0, 1, 3, 12]:
        v1 = R_all(s, m)[m]
        v2 = R_quad(s, m)
        print("s=%.5f m=%2d  R=%.12e  quad=%.12e  diff=%.2e" % (s, m, v1, v2, v1-v2))

print("=== S2: R_m(s->1) -> 0 (kernel vanishes at x=0) ===")
for s in [0.9999, 0.999999]:
    R = R_all(s, 12)
    print("s=%.7f  max|R_m(0..12)| = %.3e" % (s, np.max(np.abs(R))))

print("=== S3: operator test: f = x^2(1-x)^3 = s^6 - 2 s^8 + s^10 ===")
def ftest_s(s):
    return s**6 - 2.0*s**8 + s**10
def Aftest_exact(x):
    return (-18.0*T_moment(x, 0) + 72.0*T_moment(x, 1) - 60.0*T_moment(x, 2))
# f''' in s-coords: c_k s^{k-6}, k=6,8,10
def Af_s(s):
    R = R_all(s, 5)
    return (2.0*c_coeff(6)*R[1] - 2.0*2.0*c_coeff(8)*R[3] + 2.0*c_coeff(10)*R[5])
for s in [0.1, 0.35, 0.6, 0.8, 0.95, 0.999]:
    x = 1.0-s*s
    v1 = Af_s(s)
    v2 = Aftest_exact(x)
    print("s=%.3f x=%.4f  Af(s-method)=%.12e  exact=%.12e  diff=%.2e" % (s, x, v1, v2, v1-v2))

print("=== S4: quick eigenvalue convergence ===")
from scipy.linalg import eig as geig
for N in [10, 16, 24, 32, 44]:
    ss, Th, Fm = build_s_matrices(N)
    w, Vr = geig(Th, Fm, homogeneous_eigvals=True)
    alpha, beta = w
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = alpha/beta
    ok = np.isfinite(lam) & (np.abs(beta) > 1e-9*np.max(np.abs(alpha))) \
        & (np.abs(lam.imag) < 1e-7*np.maximum(1.0, np.abs(lam.real))) & (lam.real < 0)
    lams = np.sort(-lam[ok].real)
    a = Vr[:, ok][:, np.argsort(-lam[ok].real)][:, 0].real
    a = a/(-np.sum(a))          # f(x=1) = -sum a = 1
    fh = f_from_coeffs(a, 0.5)
    print("N=%3d lam1=%.12f lam2=%.10f lam3=%.8f f(1/2)=%.12f" % (N, lams[0], lams[1], lams[2], fh))
print("ALL SLIB TESTS DONE")
