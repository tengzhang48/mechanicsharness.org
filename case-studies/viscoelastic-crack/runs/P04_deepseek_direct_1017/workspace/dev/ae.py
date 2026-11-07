import numpy as np
from scipy import special, integrate
import time

def A(q):
    return np.exp(0.5*np.log(np.pi)+special.loggamma(1-1j*q)-special.loggamma(0.5-1j*q))
def Mminus(q): return A(-q)

def I2(K=400):
    f=lambda k: 1/(1j*k*Mminus(k))+1j/k
    # integrand ~ k^-3/2, integrate with substitution k=u^2 maybe
    def g(u):
        k=u*u
        return f(k)*2*u
    re=integrate.quad(lambda u: g(u).real,0,np.sqrt(K),limit=800,epsabs=1e-13,epsrel=1e-13)[0]
    re+=integrate.quad(lambda u: g(-u).real,0,np.sqrt(K),limit=800,epsabs=1e-13,epsrel=1e-13)[0]
    im=integrate.quad(lambda u: g(u).imag,0,np.sqrt(K),limit=800,epsabs=1e-13,epsrel=1e-13)[0]
    im+=integrate.quad(lambda u: g(-u).imag,0,np.sqrt(K),limit=800,epsabs=1e-13,epsrel=1e-13)[0]
    return re+1j*im

def I1(eps,K=400):
    # int (e^{ik eps}-1)/(ik M_-) dk = int e^{ik eps} G - int G, G=1/(ikM_-)
    # regularize G at 0: G = G0 + (G-G0), G0 =? 1/M_-(0)=1, G~1/(ik)
    # use symmetric integration and subtract
    def G(k):
        if abs(k)<1e-12: return 0j
        return 1/(1j*k*Mminus(k))
    # the -1 part is -I2? int G = I2 -? I2=int(G+ i/k)=int G + i PV int 1/k = int G (PV zero)
    # so int G = I2
    I2v=I2(K)
    # F(eps)=int e^{ik eps} G(k) dk = 2 int_0^inf Re? G(-k)=conj? 
    # use 2*Re int_0^inf e^{ik eps}G(k) for real result? check symmetry
    def h(k):
        return np.exp(1j*k*eps)*G(k)
    # integrate over [0,K], then use h(-k)=conj(h(k)) if G(-k)=conj(G(k))
    val=2*integrate.quad(lambda k: h(k).real,0,K,limit=2000,epsabs=1e-12,epsrel=1e-12)[0]
    return val - I2v

if __name__=='__main__':
    t0=time.time()
    I2v=I2()
    print('I2',I2v, 'time',time.time()-t0)
    for eps in [0.02,0.05,0.1,0.136,0.141,0.2,0.3,0.5,1.0,2.0,5.0]:
        t0=time.time()
        v=I1(eps)
        a=-v/I2v
        print('eps',eps,'I1',v,'a',a,'time',time.time()-t0)
