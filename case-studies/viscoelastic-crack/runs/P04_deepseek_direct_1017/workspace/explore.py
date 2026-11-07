import numpy as np
from scipy import special, integrate

# kernel factor M(q)=pi q coth(pi q) R(q) = A(q) R(q) A(-q)
# A(q)=sqrt(pi) Gamma(1-iq)/Gamma(1/2-iq)
# M_-(q) = A(-q) = sqrt(pi) Gamma(1+iq)/Gamma(1/2+iq)

def A(q):
    # valid complex q, use loggamma for stability
    return np.exp(0.5*np.log(np.pi) + special.loggamma(1-1j*q) - special.loggamma(0.5-1j*q))

def Mminus(q):
    return A(-q)

def R(q,V,lam):
    return (1-1j*V*q)/(1-1j*lam*V*q)

def M(q,V,lam):
    # pi q coth(pi q) R
    return np.pi*q/np.tanh(np.pi*q)*R(q,V,lam)

def integrand_I1(k, eps, V, lam):
    # (exp(i k eps)-1)/(i k M_-(k)), regular at 0
    return (np.exp(1j*k*eps)-1)/(1j*k*Mminus(k))

def integrand_I2reg(k, V, lam):
    # 1/(i k M_-) + i/k  (regular at 0, decays)
    if k == 0:
        return 0.0
    return 1/(1j*k*Mminus(k)) + 1j/k

def I1(eps, V, lam, c=0.2, K=200):
    # integrate along Im k = c (real axis shifted up); no poles of integrand in 0<Im<0.5
    f = lambda t: integrand_I1(t+1j*c, eps, V, lam)
    re = integrate.quad(lambda t: f(t).real, -K, K, limit=500, epsabs=1e-13, epsrel=1e-13)[0]
    im = integrate.quad(lambda t: f(t).imag, -K, K, limit=500, epsabs=1e-13, epsrel=1e-13)[0]
    return re+1j*im

def I2(V, lam, K=200):
    f = lambda t: integrand_I2reg(t, V, lam)
    # avoid the removable point at 0
    re = integrate.quad(lambda t: f(t).real, -K, -1e-8, limit=500, epsabs=1e-13, epsrel=1e-13)[0]
    re += integrate.quad(lambda t: f(t).real, 1e-8, K, limit=500, epsabs=1e-13, epsrel=1e-13)[0]
    im = integrate.quad(lambda t: f(t).imag, -K, -1e-8, limit=500, epsabs=1e-13, epsrel=1e-13)[0]
    im += integrate.quad(lambda t: f(t).imag, 1e-8, K, limit=500, epsabs=1e-13, epsrel=1e-13)[0]
    return re+1j*im

if __name__ == '__main__':
    lam=0.1
    for V in [0.0,0.001,0.1,10.0]:
        I2v=I2(V,lam)
        print('V',V,'I2',I2v)
        for eps in [0.1,0.5,1.0,2.0,5.0]:
            I1v=I1(eps,V,lam)
            a=-I1v/I2v
            print('  eps',eps,'I1',I1v,'a',a)
