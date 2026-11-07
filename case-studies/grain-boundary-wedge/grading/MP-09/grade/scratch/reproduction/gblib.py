"""Core library for the grain-boundary stress-relaxation problem.

Dimensionless formulation (x = z/h in (0,1), tau = t/t0):
    d sigma/d tau = (A sigma)(x) = PV int_0^1 S(x,xi) sigma_{xixixi}(xi,tau) dxi
    S(x,xi) = 1/(x-xi) - 1/(x+xi) - 2 xi (x-xi)/(x+xi)^3
BCs: sigma(0)=0, sigma'(1)=0, sigma''(1)=0.  Initial: sigma(x,0)=1 (units of sigma0).
Eigenproblem: A f = -lambda f.
"""
import numpy as np

# ---------------------------------------------------------------- kernel
def S_kernel(x, z):
    """S(x,z), singular at x=z (never evaluate there)."""
    return 1.0/(x-z) - 1.0/(x+z) - 2.0*z*(x-z)/(x+z)**3

def S_int(x, z):
    """Antiderivative of S in z (PV-safe via logs). d/dz S_int = S. Verified analytically."""
    return (-np.log(np.abs(x-z)) + np.log(x+z)
            + 6.0*x/(x+z) - 2.0*x**2/(x+z)**2)

# -------------------------------------------------- monomial moments (tests)
def T_moment(x, m):
    """T_m(x) = PV int_0^1 xi^m S(x,xi) dxi, via stable recurrences."""
    x = float(x)
    P = np.zeros(m+1); I = np.zeros(m+3); K = np.zeros(m+3); J = np.zeros(m+3)
    P[0] = -np.log((1.0-x)/x)
    for mm in range(1, m+1):
        P[mm] = -np.sum(x**np.arange(mm)/(mm-np.arange(mm))) - x**mm*np.log((1.0-x)/x)
    I[0] = np.log((1.0+x)/x)
    for mm in range(1, m+3):
        I[mm] = 1.0/mm - x*I[mm-1]
    K[0] = 1.0/x - 1.0/(1.0+x)
    for mm in range(1, m+3):
        K[mm] = I[mm-1] - x*K[mm-1]
    J[0] = 0.5*(1.0/x**2 - 1.0/(1.0+x)**2)
    for mm in range(1, m+3):
        J[mm] = K[mm-1] - x*J[mm-1]
    return P[m] - I[m] - 2.0*x*J[m+1] + 2.0*J[m+2]

def A_monomial(x, j):
    """(A x^j)(x) = j(j-1)(j-2) T_{j-3}(x); 0 for j<3."""
    if j < 3:
        return 0.0*x if np.ndim(x) else 0.0
    return j*(j-1)*(j-2)*np.array([T_moment(xx, j-3) for xx in np.atleast_1d(x)])

# ------------------------------------------- Chebyshev moment machinery
def chebT(k, u):
    """T_k(u), array-safe."""
    return np.cos(k*np.arccos(np.clip(u, -1.0, 1.0)))

def Q_int(k):
    """Q_k = int_{-1}^{1} T_k(u) du."""
    if k % 2 == 1:
        return 0.0
    return 2.0/(1.0-k*k)

def C_cauchy(v, N):
    """C_k(v) = PV int_{-1}^{1} T_k(u)/(v-u) du for k=0..N, |v|<1.
    Forward recurrence is neutrally stable on (-1,1)."""
    C = np.zeros(N+1)
    C[0] = np.log((1.0+v)/(1.0-v))
    if N >= 1:
        C[1] = v*C[0] - 2.0
    for k in range(1, N):
        C[k+1] = 2.0*v*C[k] - C[k-1] - 2.0*Q_int(k)
    return C

def D_image(vp, p, N, Kmax=None):
    """D_k^(p) = int_{-1}^{1} T_k(u)/(vp+u)^p du, vp>1, k=0..N.
    3-term relation D_{k-1} + 2 vp D_k + D_{k+1} = 2 D_k^{(p-1)}
    solved as a diagonally dominant tridiagonal system (2 vp > 2)."""
    if Kmax is None:
        rho = vp - np.sqrt(vp*vp-1.0)
        Kmax = int(N + 60 + np.log(1e-17)/np.log(rho)) + 10
    if p == 1:
        Dprev = np.array([Q_int(k) for k in range(Kmax+2)])
        D0 = np.log((vp+1.0)/(vp-1.0))
    else:
        Dprev = D_image(vp, p-1, Kmax+1, Kmax=Kmax+1)
        if p == 2:
            D0 = 2.0/(vp*vp-1.0)
        else:
            D0 = 2.0*vp/(vp*vp-1.0)**2
    K = Kmax
    M = K
    b = 2.0*vp*np.ones(M)
    c = np.ones(M-1)
    a = np.ones(M-1)
    d = 2.0*Dprev[1:K+1].copy()
    d[0] -= D0
    for i in range(1, M):
        w = a[i-1]/b[i-1]
        b[i] -= w*c[i-1]
        d[i] -= w*d[i-1]
    xvec = np.zeros(M)
    xvec[M-1] = d[M-1]/b[M-1]
    for i in range(M-2, -1, -1):
        xvec[i] = (d[i]-c[i]*xvec[i+1])/b[i]
    out = np.zeros(N+1)
    out[0] = D0
    n_take = min(N, K)
    out[1:n_take+1] = xvec[:n_take]
    return out

def _tri_solve(vp, d, D0):
    """Solve D_{k-1} + 2 vp D_k + D_{k+1} = 2 d_k (k=1..K) for D_1..D_K,
    D_0 = D0 given, D_{K+1} = 0.  d = array of D_k^{(p-1)} values, k=1..K."""
    K = len(d)
    b = 2.0*vp*np.ones(K)
    cc = np.ones(K-1)
    aa = np.ones(K-1)
    dd = 2.0*d.copy()
    dd[0] -= D0
    for i in range(1, K):
        w = aa[i-1]/b[i-1]
        b[i] -= w*cc[i-1]
        dd[i] -= w*dd[i-1]
    xv = np.zeros(K)
    xv[K-1] = dd[K-1]/b[K-1]
    for i in range(K-2, -1, -1):
        xv[i] = (dd[i]-cc[i]*xv[i+1])/b[i]
    return xv

def D_image_all(vp, N, Kmax=None):
    """Return (D1,D2,D3): D_k^(p) = int_{-1}^1 T_k(u)/(vp+u)^p du, k=0..N."""
    if Kmax is None:
        rho = vp - np.sqrt(vp*vp-1.0)
        Kmax = int(N + 60 + np.log(1e-17)/np.log(rho)) + 10
    Q = np.array([Q_int(k) for k in range(Kmax+1)])
    # p = 1
    D0_1 = np.log((vp+1.0)/(vp-1.0))
    d1 = _tri_solve(vp, Q[1:Kmax+1], D0_1)
    D1 = np.concatenate(([D0_1], d1[:N]))
    # p = 2
    D0_2 = 2.0/(vp*vp-1.0)
    rhs2 = np.concatenate(([D0_1], d1))          # D_k^1, k=0..Kmax
    d2 = _tri_solve(vp, rhs2[1:Kmax+1], D0_2)
    D2 = np.concatenate(([D0_2], d2[:N]))
    # p = 3
    D0_3 = 2.0*vp/(vp*vp-1.0)**2
    rhs3 = np.concatenate(([D0_2], d2))
    d3 = _tri_solve(vp, rhs3[1:Kmax+1], D0_3)
    D3 = np.concatenate(([D0_3], d3[:N]))
    return D1, D2, D3

def Theta_k(x, N):
    """Theta_k(x) = PV int_0^1 S(x,xi) T_k(2 xi - 1) dxi, k=0..N."""
    v = 2.0*x - 1.0
    vp = 2.0*x + 1.0
    rho = vp - np.sqrt(vp*vp-1.0)
    Kmax = int(N + 60 + np.log(1e-17)/np.log(rho)) + 10
    C = C_cauchy(v, N)
    D1, D2, D3 = D_image_all(vp, N, Kmax=Kmax)
    return C + D1 - 12.0*x*D2 + 16.0*x*x*D3

# ----------------------------------------- triple integration with BCs
def A_int0(k, v):
    """int_{-1}^{v} T_k(u) du (exact)."""
    if k == 0:
        return v + 1.0
    if k == 1:
        return 0.25*(chebT(2, v) - 1.0)
    return 0.5*(chebT(k+1, v)/(k+1) - chebT(k-1, v)/(k-1)) \
        - 0.5*((-1.0)**(k+1)/(k+1) - (-1.0)**(k-1)/(k-1))

def A_int1(k, v):
    """int_{-1}^{v} u T_k(u) du."""
    if k == 0:
        return 0.5*(v*v - 1.0)
    return 0.5*(A_int0(k+1, v) + A_int0(k-1, v))

def A_int2(k, v):
    """int_{-1}^{v} u^2 T_k(u) du."""
    if k == 0:
        return (v**3 + 1.0)/3.0
    if k == 1:
        return 0.5*(A_int1(2, v) + A_int1(0, v))
    return 0.25*(A_int0(k+2, v) + 2.0*A_int0(k, v) + A_int0(k-2, v))

def W_k(x, k):
    """W_k(x) = int_0^x (x-xi)^2/2 T_k(2xi-1) dxi = (1/16) int_{-1}^{v} (v-u)^2 T_k du."""
    v = 2.0*x - 1.0
    return (v*v*A_int0(k, v) - 2.0*v*A_int1(k, v) + A_int2(k, v))/16.0

def G0(k):
    """int_0^1 T_k(2xi-1) dxi."""
    return 0.5*A_int0(k, 1.0)

def G1(k):
    """int_0^1 xi T_k(2xi-1) dxi."""
    return 0.25*(A_int1(k, 1.0) + A_int0(k, 1.0))

def F_k(x, k):
    """F_k(x): triple integral of T_k(2xi-1) with F(0)=0, F'(1)=F''(1)=0."""
    return W_k(x, k) + x*G1(k) - 0.5*x*x*G0(k)

