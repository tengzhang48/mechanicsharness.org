"""Model reduction and validation of the singular-integral machinery.

Dimensionless problem (x=z/h, s=t/t0, t0 = 4 pi kT h^3/(E* d_gb D_gb Omega)):

  sigma_s = int_0^1 S(x,xi) sigma_xiuiui(xi) dxi   (3rd derivative in xi),
  S(x,xi) = 1/(x-xi) - 1/(x+xi) - 2 xi (x-xi)/(x+xi)^3  (PV at xi=x)
  BC: sigma(0)=0; sigma_x(1)=sigma_xx(1)=0; IC sigma(x,0)=1.

Derived in ANSWER.md from: mu=-sigma Omega, j=(dD/kT) sigma_z, 2 u_t = -Omega j_z,
sigma = sigma0 - (E*/2pi) int S u_zeta dzeta.

This script validates:
 (1) Gauss-Jacobi quadrature vs exact beta-function moments
 (2) closed form J0(x) = PV int_0^1 sqrt(1-xi)/(x-xi) dxi vs brute quadrature
 (3) closed form int_0^1 S(x,xi) dxi = ln((1+x)/(1-x)) - 4 + 6x/(1+x) - 2x^2/(1+x)^2
     (hand-derived) vs the assembled operator applied to rho=sqrt(1-xi)*q, q=1/sqrt(1-xi)
"""
import numpy as np
import scipy.special as sp


def gauss_jacobi01(M, beta):
    """nodes/weights for int_0^1 g(x) (1-x)^beta dx ~ sum_j w_j g(x_j), beta=+-1/2.

    Golub-Welsch on Jacobi weight (1-t)^beta (1+t)^0, t in [-1,1]."""
    ab = beta
    n = np.arange(M)
    an = (0.0 - ab * ab) / ((2 * n + ab) * (2 * n + ab + 2.0))
    nn = np.arange(1, M)
    bn = 4 * nn * (nn + ab) * nn * (nn + ab) / (((2 * nn + ab) ** 2) * ((2 * nn + ab) ** 2 - 1.0))
    b0 = 2.0 ** (ab + 1.0) / (ab + 1.0)  # int_{-1}^{1} (1-t)^ab dt
    J = np.diag(an)
    J += np.diag(np.sqrt(bn), 1) + np.diag(np.sqrt(bn), -1)
    tv, V = np.linalg.eigh(J)
    w = b0 * V[0] ** 2
    x = (tv + 1.0) / 2.0
    w = w * 2.0 ** (-beta - 1.0)
    return x, w


def cheb_basis(x, N):
    """T_n(2x-1), n=0..N, x array -> shape (len(x), N+1)."""
    y = 2.0 * np.asarray(x) - 1.0
    out = np.zeros((y.size, N + 1))
    out[:, 0] = 1.0
    if N >= 1:
        out[:, 1] = y
    for n in range(1, N):
        out[:, n + 1] = 2 * y * out[:, n] - out[:, n - 1]
    return out


def J0_sqrt(x):
    """PV int_0^1 sqrt(1-xi)/(x-xi) dxi = 2 + d ln((1-d)/(1+d)), d=sqrt(1-x)."""
    d = np.sqrt(1.0 - x)
    return 2.0 + d * np.log((1.0 - d) / (1.0 + d))


def J0_invsqrt(x):
    """PV int_0^1 (1-xi)^(-1/2)/(x-xi) dxi = (1/d) ln((1-d)/(1+d)), d=sqrt(1-x)."""
    d = np.sqrt(1.0 - x)
    return np.log((1.0 - d) / (1.0 + d)) / d


def kernel_op(basis_at_nodes, nodes, weights, x_col, basis_at_col, J0, beta_sign):
    """K[rho](x_i), rho(xi) = w(xi)*q_n(xi) with w=(1-xi)^beta_sign(1/2 sign),
    basis_at_nodes[j,n] = q_n(nodes_j); returns A[i,n] = int_0^1 S(x_i,xi) rho_n(xi) dxi."""
    Xc = x_col[:, None]
    Xn = nodes[None, :]
    dxi = Xc - Xn  # (i,j)
    Lx = np.log(Xc / (1.0 - Xc))  # PV int dxi/(x-xi)
    # PV part: q_n(x) J0(x) + int w(xi) (q_n(xi)-q_n(x))/(x-xi) dxi
    A = (basis_at_col * J0[:, None])
    num = basis_at_nodes[None, :, :] - basis_at_col[:, None, :]  # (i,j,n)
    A = A + np.einsum('j,ijn->in', weights, num / dxi[:, :, None])
    # image parts (regular for x>0)
    k1 = -1.0 / (Xc + Xn)
    k2 = -2.0 * Xn * dxi / (Xc + Xn) ** 3
    A = A + np.einsum('j,ij,jn->in', weights, k1, basis_at_nodes)
    A = A + np.einsum('j,ij,jn->in', weights, k2, basis_at_nodes)
    return A


# ---------------- checks ----------------
rng = np.random.default_rng(0)
print("=== check 1: Gauss-Jacobi moments vs Beta function ===")
M = 120
xg, wg = gauss_jacobi01(M, 0.5)
for k in range(0, 6):
    mine = np.sum(wg * xg ** k)
    exact = sp.beta(k + 1.0, 1.5)
    print(f"  k={k}: quad={mine:.16e} beta={exact:.16e} rel={abs(mine-exact)/exact:.2e}")
xg2, wg2 = gauss_jacobi01(M, -0.5)
for k in range(0, 6):
    mine = np.sum(wg2 * xg2 ** k)
    exact = sp.beta(k + 1.0, 0.5)
    print(f"  invsqrt k={k}: quad={mine:.16e} beta={exact:.16e} rel={abs(mine-exact)/exact:.2e}")

print("=== check 2: J0 closed form vs brute PV quadrature ===")


def brute_pv(x, g):
    """PV int_0^1 g(xi)/(x-xi) dxi by subtraction with fine Gauss-Legendre."""
    t, w = np.polynomial.legendre.leggauss(2000)
    xi = (t + 1) / 2
    ww = w / 2
    gs = g(xi)
    gx = g(x)
    integrand = np.interp(xi, xi, (gs - gx)) / (x - xi)
    return np.sum(ww * np.interp(xi, xi, integrand)) + gx * np.log(x / (1 - x))


for x in [0.05, 0.2, 0.5, 0.8, 0.95]:
    g = lambda v: np.sqrt(1 - v)
    exact = J0_sqrt(x)
    b = brute_pv(x, g)
    print(f"  x={x}: J0={exact:.14f} brute={b:.14f} diff={exact-b:.2e}")
for x in [0.05, 0.2, 0.5, 0.8, 0.95]:
    g = lambda v: 1.0 / np.sqrt(1 - v)
    exact = J0_invsqrt(x)
    b = brute_pv(x, g)
    print(f"  invsqrt x={x}: J0={exact:.14f} brute={b:.14f} diff={exact-b:.2e}")

print("=== check 3: int_0^1 S(x,xi) dxi vs hand-derived closed form ===")
# rho(xi) = 1 = sqrt(1-xi) * q(xi) with q = (1-xi)^(-1/2): not polynomial...
# instead check with rho(xi)=sqrt(1-xi)*q, q=1/sqrt(1-xi)?? Simpler: check
# int S(x,xi) dxi directly by brute quadrature and compare to closed form.
def closed_intS(x):
    return (np.log((1 + x) / (1 - x)) - 4 + 6 * x / (1 + x) - 2 * x ** 2 / (1 + x) ** 2)

t, w = np.polynomial.legendre.leggauss(4000)
xi = (t + 1) / 2
ww = w / 2
for x in [0.01, 0.05, 0.1, 0.3, 0.6, 0.9, 0.99]:
    S = 1.0 / (x - xi) - 1.0 / (x + xi) - 2 * xi * (x - xi) / (x + xi) ** 3
    pv = np.log(x / (1 - x))
    total = np.sum(ww * S) - np.sum(ww) / (x - x) * 0  # naive: replace PV part
    # do PV by subtraction only for the 1/(x-xi) piece
    naive = np.sum(ww * S)
    sub = np.sum(ww * (1.0 - 1.0) / (x - xi))
    brute = np.sum(ww * (S - 1.0 / (x - xi))) + pv
    print(f"  x={x}: brute={brute:.14f} closed={closed_intS(x):.14f} diff={brute-closed_intS(x):.2e}")

print("=== check 2b: J0_invsqrt closed form vs t-substituted brute ===")
t, w = np.polynomial.legendre.leggauss(4000)
for x in [0.05, 0.2, 0.5, 0.8, 0.95]:
    d2 = 1.0 - x
    # PV 2*int_0^1 dt/(t^2-d^2): subtract value at t^2=d^2 -> integrand (t^2-d^2)/(t^2-d^2)=1
    tt = (t + 1) / 2
    ww = w / 2
    f = 2.0 / (tt ** 2 - d2)
    g = tt ** 2 - d2  # f*g = 2
    # PV: 2/(t^2-d^2) = 2 (t^2-d^2)/(t^2-d^2)^2 -> split: 1 + d2/(t^2-d^2)...
    # do: PV int = int [2/(t^2-d^2) - 2/(2 t d?)]. Simplest exact: analytic.
    exact = J0_invsqrt(x)
    # numeric PV via splitting the interval at t=sqrt(d2) and subtraction of 1/(t-d)-type singularity
    d = np.sqrt(d2)
    def h(u):  # u in place of t, with singular 1/(u-d)
        return 2.0 / (u ** 2 - d ** 2) - 2.0 / (2 * d * (u - d))
    def hs(u):
        return 1.0 / (u + d)
    # 2/(t^2-d^2) = 1/(d(t-d)) - 1/(d(t+d))
    def fa(u):
        return -1.0 / (d * (u + d))
    from numpy.polynomial.legendre import leggauss
    def quad(lo, hi, func, n=2000):
        z, wz = leggauss(n)
        m = (lo + hi) / 2; r = (hi - lo) / 2
        return r * np.sum(wz * func(m + r * z))
    pv = quad(0, d, fa) + quad(d, 1, fa)  # the -1/(d(t+d)) part is regular
    pv += (1.0 / d) * (np.log(d) - np.log(0) * 0)  # PV int_0^? ... handle below
    # PV int_0^1 dt/(t-d) = ln(1-d) - ln(d)
    pv += (1.0 / d) * (np.log(1 - d) - np.log(d))
    print(f"  x={x}: J0_invsqrt={exact:.14f} brute={pv:.14f} diff={exact-pv:.2e}")
