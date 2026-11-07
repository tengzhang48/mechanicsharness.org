import numpy as np
from scipy import special, integrate
from whfft import A,Mm,Mp,S0

def res1(n):
    return (-1)**n/(np.sqrt(np.pi)*special.factorial(n)*1j*special.gamma(0.5-n))

def Gplus_resum(q,eps,N=60):
    # G_+ = G(q) - eps/Mm(q) - sum res/(s(kn-q)) + sum res e^{-s eps}/(s(kn-q))
    G=S0(q,eps)/Mm(q)
    s=G - eps/Mm(q)
    for n in range(N):
        sn=n+0.5; kn=1j*sn; r=res1(n)
        s += r*(np.exp(-sn*eps)-1)/(sn*(kn-q))
    return s

def Gplus_direct(q,eps,K=400):
    def G(k): return S0(k,eps)/Mm(k)
    f=lambda k: G(k)/(k-q)
    val=0.5*G(q)
    re=integrate.quad(lambda k: f(k).real,-K,q-1e-6,limit=20000,epsabs=1e-13,epsrel=1e-12)[0]
    re+=integrate.quad(lambda k: f(k).real,q+1e-6,K,limit=20000,epsabs=1e-13,epsrel=1e-12)[0]
    im=integrate.quad(lambda k: f(k).imag,-K,q-1e-6,limit=20000,epsabs=1e-13,epsrel=1e-12)[0]
    im+=integrate.quad(lambda k: f(k).imag,q+1e-6,K,limit=20000,epsabs=1e-13,epsrel=1e-12)[0]
    val+=(re+1j*im)/(2*np.pi*1j)
    return val

if __name__=='__main__':
    # partial fraction identity
    for q in [0.7,1.5,3.0]:
        sm=sum(res1(n)/(1j*(n+0.5)-q) for n in range(2000))
        print('pf',q,sm,'-1/Mm',-1/Mm(q))
    eps=0.125
    for q in [0.5,1.5,3.0]:
        print('q',q)
        for N in [10,20,40,80]:
            print('  resum N',N,Gplus_resum(q,eps,N))
        print('  direct',Gplus_direct(q,eps,K=400))
