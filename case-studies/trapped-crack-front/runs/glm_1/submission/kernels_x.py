"""
X-domain 2.5D traction kernels (closed form) and finite-part application to DD basis
functions on the crack X>0 (units 1/|k|, mode k=+1).

T[b](X') = int_0^inf K(X'-X) b(X) dX (Hadamard finite part at X'=X),
K_ij(t) = (1/2pi) int M_ij(p,1) e^{ipt} dp  (mu=1):
  K_yy(t) =  [K1(|t|)/|t|] /(2 pi (1-nu))
  K_xx(t) =  [K1(|t|)/|t| + nu K0(|t|)] /(2 pi (1-nu))
  K_zz(t) =  [K1(|t|)/|t| - nu K0(|t|)/(1-nu)] /(2 pi)
  K_xz(t) = -i nu sgn(t) K1(|t|) /(2 pi (1-nu))
Small-|t|: K1(|t|)/|t| = 1/t^2 + (1/2)[ln(|t|/2)+gamma-1/2] + O(t^2 ln t),
           K0(|t|)     = -ln(|t|/2) - gamma + O(t^2).
On-crack evaluation subtracts the local model phi = exp(-(X-Xj)):
  I1 = PV int b/(Xj-X)    = int (b - b0 phi)/(Xj-X) + b0 Ei(Xj)
  I2 = fp  int b/(Xj-X)^2 = int (b - S)/(Xj-X)^2 + b0 e^{Xj}(J-1/Xj) - b1 e^{Xj} J + b2/2
       S = phi (b0 + b1 (X-Xj) + b2 (X-Xj)^2 / 2),  J = e^{-Xj} Ei(Xj)
(model integrals: A2 = e^{Xj}(J - 1/Xj), A3 = -e^{Xj} J, A4 = 1.)
Ahead evaluation (X'<0) is nonsingular (validated against contour quadrature).
"""
import numpy as np
from scipy.special import k0 as besselk0, k1 as besselk1, expi

GAMMA = 0.5772156649015328606
NU = 0.3
MU = 1.0

def kernels(t):
    at = np.abs(np.asarray(t, float)); sgn = np.sign(np.asarray(t, float))
    K1ot = besselk1(at)/at
    K0 = besselk0(at)
    Kyy = MU/(2*np.pi*(1-NU))*K1ot
    Kxx = MU/(2*np.pi*(1-NU))*(K1ot + NU*K0)
    Kzz = MU/(2*np.pi)*(K1ot - NU/(1-NU)*K0)
    Kxz = -1j*MU*NU/(2*np.pi*(1-NU))*sgn*besselk1(at)
    return Kyy, Kxx, Kzz, Kxz

def singular_coeffs():
    c_half = 0.5*(GAMMA - 0.5 - np.log(2.0))
    c_k0 = np.log(2.0) - GAMMA
    A_yy = MU/(2*np.pi*(1-NU)); B_yy = MU/(4*np.pi*(1-NU)); C_yy = MU/(2*np.pi*(1-NU))*c_half
    A_xx = MU/(2*np.pi*(1-NU)); B_xx = MU/(2*np.pi*(1-NU))*(0.5 - NU)
    C_xx = MU/(2*np.pi*(1-NU))*(c_half + NU*c_k0)
    A_zz = MU/(2*np.pi); B_zz = MU/(4*np.pi) + MU*NU/(2*np.pi*(1-NU))
    C_zz = MU/(2*np.pi)*(c_half - NU/(1-NU)*c_k0)
    D_xz = -1j*MU*NU/(2*np.pi*(1-NU))
    return {'yy': (A_yy, 0.0, B_yy, C_yy), 'xx': (A_xx, 0.0, B_xx, C_xx),
            'zz': (A_zz, 0.0, B_zz, C_zz), 'xz': (0.0, D_xz, 0.0, 0.0)}

def K_analytic(t, entry):
    A, D, B, C = singular_coeffs()[entry]
    at = np.abs(np.asarray(t, float)); t = np.asarray(t, float)
    Kyy, Kxx, Kzz, Kxz = kernels(t)
    K = {'yy': Kyy, 'xx': Kxx, 'zz': Kzz, 'xz': Kxz}[entry]
    return K - (A/at**2 + B*np.log(at) + C + D/t)

class BasisFun:
    """b(X) = e^{-X} sum_k c_k X^{alpha+k} (Laguerre), divided by `scale`.
    Edge behavior: b ~ X^alpha / scale as X->0 (L_m(0)=1)."""
    def __init__(self, alpha, m, scale=1.0):
        from math import comb, factorial
        c = np.zeros(m+1)
        for i in range(m+1):
            c[i] = (-1.0)**i * comb(m, i) / factorial(i)
        self.alpha = alpha; self.c = c/scale; self.m = m; self.scale = scale
    def eval(self, X):
        X = np.asarray(X, float); k = np.arange(len(self.c))
        return np.exp(-X)*(np.power.outer(X, self.alpha+k) @ self.c)
    def derivs(self, X):
        k = np.arange(len(self.c)); c = self.c; a = self.alpha
        e = np.exp(-X)
        pw = X**(a+k); pw1 = X**(a+k-1); pw2 = X**(a+k-2)
        b  = e*(pw @ c)
        d1 = e*((a+k)*pw1 @ c) - b
        d2 = e*(((a+k)*(a+k-1)*pw2) @ c) - 2*e*((a+k)*pw1 @ c) + b
        return b, d1, d2
    def moment(self):
        from scipy.special import gamma as G
        k = np.arange(len(self.c))
        return float((G(self.alpha+k+1)*self.c).sum())

def Jfun(X):
    return np.exp(-X)*expi(X)

# ---------- fixed-panel Gauss-Legendre ----------
_GLN = 20
_xg, _wg = np.polynomial.legendre.leggauss(_GLN)

def _panels(a, b, n, cluster='cubic'):
    if cluster == 'cubic':
        return a + (b-a)*(np.linspace(0, 1, n+1))**3
    if cluster == 'geom':
        return a*(b/a)**(np.linspace(0, 1, n+1))
    return np.linspace(a, b, n+1)

def nodes_for(Xj, nL=10, nR=10, Xmax=60.0):
    eL = _panels(0.0, Xj, nL, 'cubic')
    eR = _panels(Xj, Xmax, nR, 'geom')
    nodes = []; wts = []
    for e in zip(eL[:-1], eL[1:]):
        nodes.append(0.5*(e[1]+e[0]) + 0.5*(e[1]-e[0])*_xg)
        wts.append(0.5*(e[1]-e[0])*_wg)
    for e in zip(eR[:-1], eR[1:]):
        nodes.append(0.5*(e[1]+e[0]) + 0.5*(e[1]-e[0])*_xg)
        wts.append(0.5*(e[1]-e[0])*_wg)
    return np.concatenate(nodes), np.concatenate(wts)

def apply_kernels_matrix(basis, Xj, entries):
    """T[bf](Xj) for all basis functions, entries; finite parts with the
    phi = e^{-(X-Xj)} subtraction model (see module docstring)."""
    X, W = nodes_for(Xj)
    Bmat = np.array([bf.eval(X) for bf in basis])
    b0 = np.zeros(len(basis)); d1v = np.zeros(len(basis)); d2v = np.zeros(len(basis))
    for i, bf in enumerate(basis):
        b0[i], d1v[i], d2v[i] = bf.derivs(Xj)
    I0 = np.array([bf.moment() for bf in basis])
    t = Xj - X
    phi = np.exp(-(X - Xj))
    f1 = (Bmat - b0[:, None]*phi[None, :])/t[None, :]
    I1 = f1 @ W + b0*expi(Xj)
    # I2 with Hermite-matched models: phi0=e^{-u}(1+u+u^2/2), phi1=e^{-u}(u+u^2),
    # phi2=e^{-u}u^2/2 (u=X-Xj); value/slope/curvature matched at Xj, decay at inf.
    u = X - Xj
    ph0 = np.exp(-u)*(1 + u + 0.5*u*u)
    ph1 = np.exp(-u)*(u + u*u)
    ph2 = np.exp(-u)*0.5*u*u
    S = (b0[:, None]*ph0[None, :] + d1v[:, None]*ph1[None, :] + d2v[:, None]*ph2[None, :])
    f2 = (Bmat - S)/(t[None, :])**2
    ex = np.exp(Xj); Ei = expi(Xj)
    B0 = ex*(0.5 - 1.0/Xj)
    B1 = -Ei + ex
    B2 = 0.5*ex
    I2 = f2 @ W + b0*B0 + d1v*B1 + d2v*B2
    fl = Bmat*np.log(np.abs(t))[None, :]
    Ilog = fl @ W
    out = {}
    for entry in entries:
        A, D, B, C = singular_coeffs()[entry]
        Kan = K_analytic(t, entry)
        Ian = (Bmat*Kan[None, :]) @ W
        out[entry] = (A*I2 + D*I1 + B*Ilog + C*I0 + Ian)
    return out

def apply_ahead_matrix(basis, x, entries, Xmax=80.0):
    """T[bf](x) for x<0 (ahead), plain integrals (validated).
    Panels cluster down to |x|/100 to resolve the near-singular contribution."""
    xa = abs(x)
    edges = np.concatenate([xa*np.array([1e-3, 1e-2, 1e-1]),
                            xa*np.geomspace(1.0, 1e2, 5)[1:],
                            np.geomspace(100*xa, Xmax, 8)[1:]])
    if edges[-1] < Xmax:
        edges = np.concatenate([edges, [Xmax]])
    nodes = []; wts = []
    for a, b in zip(edges[:-1], edges[1:]):
        nodes.append(0.5*(a+b) + 0.5*(b-a)*_xg); wts.append(0.5*(b-a)*_wg)
    X = np.concatenate(nodes); W = np.concatenate(wts)
    Bmat = np.array([bf.eval(X) for bf in basis])
    out = {}
    for entry in entries:
        K = kernels(x - X)
        Kd = {'yy': K[0], 'xx': K[1], 'zz': K[2], 'xz': K[3]}[entry]
        out[entry] = (Bmat*Kd[None, :]) @ W
    return out

if __name__ == "__main__":
    bf = BasisFun(0.0, 0)   # b = e^{-X}
    for Xj in [0.2, 0.5, 1.0, 3.0]:
        J = Jfun(Xj)
        d = apply_kernels_matrix([bf], Xj, ['yy'])
        print(f"Xj={Xj}: (validation of I1/I2 in check_I2.py)")
