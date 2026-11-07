"""Finite crack |x|<1 in plane-strain Biot poroelasticity, Laplace domain (parameter s).
Units: L = 1, G = 1, c = 1.  Galerkin method with Chebyshev bases:
  U(x) = u_y(x,0+) = sum a_n sqrt(1-x^2) U_n(x)       (half opening)
  Q(x) = q_y(x,0+) = sum b_n T_n(x)/sqrt(1-x^2)       (flux from crack into upper solid)
Half-plane symbols (k>0, m = sqrt(k^2+s)):
  sigma_yy^ = A11 U^ + A12 Q^ ,  p^ = A21 U^ + A22 Q^
  A11 = -k/(1-nu_u) * (1 - gam*2k^2/(m(m+k))),  gam = (nu_u-nu)/(1-nu)
  A12 = -(2 beta/(1-nu_u)) k/(m(m+k)),  A21 = -s A12,  A22 = 1/(kappa m)
Loading 1: sigma_yy = -1/s, p = 0 on |x|<1.  Loading 2: sigma_yy = 0, p = 1/s.
Outputs per loading: Laplace transforms Kt (K/sqrt(pi L)), Vt (crack volume, units L^2/G).
"""
import numpy as np
from scipy.special import jv, iv, kv

class Params:
    def __init__(self, nu, nuu, alpha):
        self.nu, self.nuu, self.alpha = nu, nuu, alpha
        self.gam = (nuu - nu)/(1 - nu)
        self.beta = (nuu - nu)/(alpha*(1 - 2*nu))
        self.eta = alpha*(1 - 2*nu)/(2*(1 - nu))
        self.S = alpha**2*(1 - 2*nu)**2*(1 - nuu)/(2*(nuu - nu)*(1 - nu))  # G=1
        self.kappa = self.S  # c = 1

def gl_panels(kmax, width=0.5, npts=12):
    xg, wg = np.polynomial.legendre.leggauss(npts)
    edges = np.arange(0.0, kmax + width, width)
    a, b = edges[:-1], edges[1:]
    k = (0.5*(b - a)[:, None]*xg[None, :] + 0.5*(b + a)[:, None]).ravel()
    w = (0.5*(b - a)[:, None]*wg[None, :]).ravel()
    return k, w

def assemble(s, N, P, kmax=None, width=0.5, npts=12, ntheta=400):
    nu, nuu = P.nu, P.nuu
    n = 2*np.arange(N)                      # even Chebyshev indices
    sa = np.sqrt(complex(s)) if np.iscomplexobj(s) or (np.real(s) < 0) else np.sqrt(s)
    if kmax is None:
        kmax = max(80.0*abs(sa), 4.0*n[-1] + 200.0)
    k, w = gl_panels(kmax, width, npts)
    m = np.sqrt(k**2 + s)
    orders = np.arange(0, 2*N + 1)
    J = jv(orders[:, None], k[None, :])     # (2N+1, Nq)
    Jn = J[n, :]; Jn1 = J[n + 1, :]
    # I11 = int A11 J_{m+1} J_{n+1} / k^2 ; A11 = A11inf + R11, A11inf = -k/(1-nu)
    R11 = -(P.gam/(1 - nuu))*k*(m - k)*(m + 2*k)/(m*(m + k))
    I11 = (Jn1*(w*R11/k**2)) @ Jn1.T
    I11 = I11 + np.diag(-(1.0/(1 - nu))/(2.0*(n + 1)))
    # I12 = int A12 J_{m+1} J_n / k
    A12 = -(2*P.beta/(1 - nuu))*k/(m*(m + k))
    I12 = (Jn1*(w*A12/k)) @ Jn.T
    I21 = -s*I12.T                           # int A21 J_m J_{n+1}/k
    # I22 = int J_m J_n/(kappa m) = (1/kappa)[int J_mJ_n/k - s int J_mJ_n/(m k (m+k))]
    R22 = -s/(m*k*(m + k))
    I22 = (Jn*(w*R22)) @ Jn.T
    base = np.zeros((N, N)); base[np.arange(1, N), np.arange(1, N)] = 1.0/(2.0*n[1:])
    I22 = (I22 + base)/P.kappa
    # (0,0) entry exactly: (2/(pi kappa)) int_0^{pi/2} I0(a cos th) K0(a cos th) dth
    u, wu = np.polynomial.legendre.leggauss(ntheta)
    u = 0.5*(u + 1); wu = 0.5*wu
    phi = 0.5*np.pi*u**2; dphi = np.pi*u*wu             # theta = pi/2 - phi
    xarg = sa*np.sin(phi)
    I00 = (2/(np.pi*P.kappa))*np.sum(dphi*iv(0, xarg)*kv(0, xarg))
    I22[0, 0] = I00
    sg = (-1.0)**((n[:, None] + n[None, :])//2)
    M11 = np.pi*(n[:, None] + 1)*(n[None, :] + 1)*sg*I11
    M12 = np.pi*(n[:, None] + 1)*sg*I12
    M21 = np.pi*(n[None, :] + 1)*sg*I21
    M22 = np.pi*sg*I22
    return np.block([[M11, M12], [M21, M22]]), n

def solve(s, N, P, **kw):
    M, n = assemble(s, N, P, **kw)
    out = {}
    for name, (sc, pc) in {'L1': (-1.0/s, 0.0), 'L2': (0.0, 1.0/s)}.items():
        rhs = np.zeros(2*N, dtype=complex if np.iscomplexobj(M) else float)
        rhs[0] = sc*np.pi/2; rhs[N] = pc*np.pi
        x = np.linalg.solve(M, rhs)
        a, b = x[:N], x[N:]
        Kt = np.sum(a*(n + 1))/(1 - P.nu)     # K/sqrt(pi L), G = 1
        Vt = np.pi*a[0]
        out[name] = dict(K=Kt, V=Vt, a=a, b=b)
    return out

if __name__ == '__main__':
    P = Params(0.2, 0.4, 0.89)
    for s in [1e-4, 1e-2, 1.0, 1e2, 1e3]:
        for N in [10, 20, 40]:
            o = solve(s, N, P)
            print(f"s={s:8.1e} N={N:3d}  sK1={s*o['L1']['K']:.6f} sV1={s*o['L1']['V']:.6f}  sK2={s*o['L2']['K']:.6f} sV2={s*o['L2']['V']:.6f}")
