import numpy as np
from scipy import special, integrate

# ---------- factor functions ----------
def logA(q):
    return 0.5*np.log(np.pi) + special.loggamma(1-1j*q) - special.loggamma(0.5-1j*q)
def A(q):
    return np.exp(logA(q))
def Mminus(q):   # A(-q)
    return A(-q)
def Rfac(q,V,lam):
    return (1-1j*V*q)/(1-1j*lam*V*q)
def Mplus(q,V,lam):
    return A(q)*Rfac(q,V,lam)

# ---------- SIF / boundedness integral giving a(eps) ----------
def integrand_I1(k, eps):
    return (np.exp(1j*k*eps)-1)/(1j*k*Mminus(k))
def integrand_I2reg(k):
    if k == 0: return 0.0
    return 1/(1j*k*Mminus(k)) + 1j/k
def I1(eps, c=0.2, K=200):
    f = lambda t: integrand_I1(t+1j*c, eps)
    re = integrate.quad(lambda t: f(t).real, -K, K, limit=800, epsabs=1e-12, epsrel=1e-12)[0]
    im = integrate.quad(lambda t: f(t).imag, -K, K, limit=800, epsabs=1e-12, epsrel=1e-12)[0]
    return re+1j*im
def I2(K=200):
    f = lambda t: integrand_I2reg(t)
    re = integrate.quad(lambda t: f(t).real, -K, -1e-8, limit=800, epsabs=1e-12, epsrel=1e-12)[0]
    re += integrate.quad(lambda t: f(t).real, 1e-8, K, limit=800, epsabs=1e-12, epsrel=1e-12)[0]
    im = integrate.quad(lambda t: f(t).imag, -K, -1e-8, limit=800, epsabs=1e-12, epsrel=1e-12)[0]
    im += integrate.quad(lambda t: f(t).imag, 1e-8, K, limit=800, epsabs=1e-12, epsrel=1e-12)[0]
    return re+1j*im

# ---------- residue machinery ----------
def dA(q, h=1e-5):
    return (A(q+h)-A(q-h))/(2*h)
def dMminus(q, h=1e-5):
    return (Mminus(q+h)-Mminus(q-h))/(2*h)

def Shat(q, eps, a):
    # ((exp(i q eps)-1)+a)/(i q)
    return ((np.exp(1j*q*eps)-1)+a)/(1j*q)

def coefficients(eps, a, V, lam, N=60):
    an = np.arange(N)+0.5
    # UHP poles of H: q=+i an
    qp = 1j*an
    Shatp = np.array([Shat(q,eps,a) for q in qp])
    dM = np.array([dMminus(q) for q in qp])
    c = Shatp/dM
    return an, qp, c

def U_open(X, eps, a, V, lam, N=60, Nsum=None):
    """Compute U(X) for X<0 via residue sum."""
    an, qp, c = coefficients(eps,a,V,lam,N)
    if Nsum is None: Nsum=N
    # LHP zeros of M_+ :
    qzero = [-1j*(n+0.5) for n in range(N)]
    if V>0:
        qzero.append(-1j/V)
    # H_reg(q)=H(q)+i a/q
    def Hreg(q):
        return Shat(q,eps,a)/Mminus(q) + 1j*a/q
    # H_+(q)=Hreg(q)-sum_n c_n/(q - i an)
    def Hplus(q):
        s = Hreg(q)
        s -= np.sum(c/(q-qp))
        return s
    # derivative M_+'
    def dMplus(q):
        h=1e-5
        return (Mplus(q+h,V,lam)-Mplus(q-h,V,lam))/(2*h)
    val = a*1.0  # constant far field
    for qj in qzero:
        val += 1j*Hplus(qj)/dMplus(qj)*np.exp(1j*qj*X)
    return val

if __name__=='__main__':
    lam=0.1; dbar=0.2
    I2v=I2()
    for V in [0.0]:
        for eps in [0.05,0.1,0.141,0.2,0.3,0.5,1.0]:
            a=-I1(eps)/I2v
            u=U_open(-eps,eps,a,V,lam,N=80)
            print('V',V,'eps',eps,'a',a,'U(-eps)',u,'need',dbar/2)
