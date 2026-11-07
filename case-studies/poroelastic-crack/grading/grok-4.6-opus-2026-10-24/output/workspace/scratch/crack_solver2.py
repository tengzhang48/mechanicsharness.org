"""Finite crack |x|<1, plane-strain Biot poroelasticity, Laplace domain (s real or complex).
Same formulation as crack_solver.py, but the k-integrals are split into an exact
Gauss-Legendre part on [0,kmax] and a tail on [kmax,inf) in which the Bessel product
J_mu J_nu is replaced by its non-oscillatory Debye form
  (1/pi) (k^2-mu^2)^(-1/4) (k^2-nu^2)^(-1/4) cos(psi_mu - psi_nu),
  psi_mu = sqrt(k^2-mu^2) - mu*arccos(mu/k).
Units L = G = c = 1.
"""
import numpy as np
from scipy.special import jv, iv, kv, ive, kve

class Params:
    def __init__(self, nu, nuu, alpha):
        self.nu, self.nuu, self.alpha = nu, nuu, alpha
        self.gam = (nuu - nu)/(1 - nu)
        self.beta = (nuu - nu)/(alpha*(1 - 2*nu))
        self.eta = alpha*(1 - 2*nu)/(2*(1 - nu))
        self.S = alpha**2*(1 - 2*nu)**2*(1 - nuu)/(2*(nuu - nu)*(1 - nu))
        self.kappa = self.S

def gl_panels(a0, kmax, width, npts, scale=None):
    """Gauss-Legendre panels on [0,kmax].  If scale (=sqrt|s|) < 1, geometric panels
    (10 per decade) resolve k ~ scale on [1e-7*scale, 1]; uniform panels beyond, of width
    0.5 if scale < 4 else `width`."""
    xg, wg = np.polynomial.legendre.leggauss(npts)
    if scale is not None and scale < 1.0:
        g = np.geomspace(1e-7*scale, 1.0, int(round(10*np.log10(1.0/(1e-7*scale)))) + 1)
        head = [0.0] + list(g[:-1]); start = 1.0
    else:
        head = []; start = 0.0
    w = 0.5 if (scale is not None and scale < 4.0) else width
    nP = max(1, int(np.ceil((kmax - start)/w)))
    edges = np.array(head + list(np.linspace(start, kmax, nP + 1)))
    a, b = edges[:-1], edges[1:]
    k = (0.5*(b - a)[:, None]*xg[None, :] + 0.5*(b + a)[:, None]).ravel()
    w = (0.5*(b - a)[:, None]*wg[None, :]).ravel()
    return k, w

def symbols(k, s, P):
    m = np.sqrt(k**2 + s)
    R11 = -(P.gam/(1 - P.nuu))*k*(m - k)*(m + 2*k)/(m*(m + k))
    A12 = -(2*P.beta/(1 - P.nuu))*k/(m*(m + k))
    R22 = -s/(m*k*(m + k))
    return R11, A12, R22

def debye_products(k, mus, nus):
    """non-oscillatory part of J_mu(k) J_nu(k), arrays (len(mus), len(nus), len(k))"""
    mu = np.asarray(mus, float)[:, None]; nu = np.asarray(nus, float)[:, None]
    sm = np.sqrt(k[None, :]**2 - mu**2); sn = np.sqrt(k[None, :]**2 - nu**2)
    pm = sm - mu*np.arccos(mu/k[None, :]); pn = sn - nu*np.arccos(nu/k[None, :])
    amp_m = sm**-0.5; amp_n = sn**-0.5
    # cos(pm - pn) = cos pm cos pn + sin pm sin pn  -> separable sums
    Cm, Sm_ = amp_m*np.cos(pm), amp_m*np.sin(pm)
    Cn, Sn_ = amp_n*np.cos(pn), amp_n*np.sin(pn)
    return Cm, Sm_, Cn, Sn_

def tail_matrix(F, k, w, mus, nus):
    Cm, Sm_, Cn, Sn_ = debye_products(k, mus, nus)
    return ((Cm*(w*F)) @ Cn.T + (Sm_*(w*F)) @ Sn_.T)/np.pi

def assemble(s, N, P, kmax=None, width=2.0, npts=16, ntheta=600, ntail=200):
    nu, nuu = P.nu, P.nuu
    n = 2*np.arange(N)
    sa = np.sqrt(s + 0j) if np.iscomplexobj(s) else np.sqrt(s)
    if kmax is None:
        kmax = 4.0*(2*N + 1) + 40.0*abs(sa) + 200.0
    k_all, w_all = gl_panels(0.0, kmax, width, npts, scale=abs(sa))
    dt = complex if np.iscomplexobj(s) else float
    I11 = np.zeros((N, N), dt); I12 = np.zeros((N, N), dt); I22 = np.zeros((N, N), dt)
    chunk = max(1000, int(4e7//(2*N + 1)))
    for i0 in range(0, k_all.size, chunk):
        k = k_all[i0:i0 + chunk]; w = w_all[i0:i0 + chunk]
        Jn = jv(n[:, None], k[None, :]); Jn1 = jv(n[:, None] + 1, k[None, :])
        R11, A12, R22 = symbols(k, s, P)
        I11 += (Jn1*(w*R11/k**2)) @ Jn1.T
        I12 += (Jn1*(w*A12/k)) @ Jn.T
        I22 += (Jn*(w*R22)) @ Jn.T
    # tails on [kmax, inf): k = kmax/u, u in (0,1]
    ug, wg = np.polynomial.legendre.leggauss(ntail)
    u = 0.5*(ug + 1); wu = 0.5*wg
    kt = kmax/u**2; wt = wu*2*kmax/u**3          # k = kmax/u^2 for faster tail sampling
    T11, T12, T22 = symbols(kt, s, P)
    I11 = I11 + tail_matrix(T11/kt**2, kt, wt, n + 1, n + 1)
    I12 = I12 + tail_matrix(T12/kt, kt, wt, n + 1, n)
    I22 = I22 + tail_matrix(T22, kt, wt, n, n)
    I11 = I11 + np.diag(-(1.0/(1 - nu))/(2.0*(n + 1)))
    I21 = -s*I12.T
    base = np.zeros((N, N)); base[np.arange(1, N), np.arange(1, N)] = 1.0/(2.0*n[1:])
    I22 = (I22 + base)/P.kappa
    ug, wg = np.polynomial.legendre.leggauss(ntheta)
    u = 0.5*(ug + 1); wu = 0.5*wg
    phi = 0.5*np.pi*u**2; dphi = np.pi*u*wu
    xarg = sa*np.sin(phi)
    # I0 K0 via exponentially scaled functions: iv*kv = ive*kve*exp(-i Im z) for Re z > 0
    prod = ive(0, xarg)*kve(0, xarg)*np.exp(-1j*np.imag(xarg)) if np.iscomplexobj(xarg) else ive(0, xarg)*kve(0, xarg)
    I22[0, 0] = (2/(np.pi*P.kappa))*np.sum(dphi*prod)
    sg = (-1.0)**((n[:, None] + n[None, :])//2)
    M11 = np.pi*(n[:, None] + 1)*(n[None, :] + 1)*sg*I11
    M12 = np.pi*(n[:, None] + 1)*sg*I12
    M21 = np.pi*(n[None, :] + 1)*sg*I21
    M22 = np.pi*sg*I22
    return np.block([[M11, M12], [M21, M22]]), n

def solve(s, N, P, **kw):
    M, n = assemble(s, N, P, **kw)
    out = {}
    dt = complex if np.iscomplexobj(M) else float
    for name, (sc, pc) in {'L1': (-1.0/s, 0.0), 'L2': (0.0, 1.0/s)}.items():
        rhs = np.zeros(2*N, dtype=dt)
        rhs[0] = sc*np.pi/2; rhs[N] = pc*np.pi
        x = np.linalg.solve(M, rhs)
        a, b = x[:N], x[N:]
        out[name] = dict(K=np.sum(a*(n + 1))/(1 - P.nu), V=np.pi*a[0], a=a, b=b)
    return out
