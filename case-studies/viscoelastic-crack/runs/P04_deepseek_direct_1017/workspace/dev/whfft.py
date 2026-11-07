import numpy as np
from scipy import special, integrate

def A(q): return np.exp(0.5*np.log(np.pi)+special.loggamma(1-1j*q)-special.loggamma(0.5-1j*q))
def Mm(q): return A(-q)
def Mp(q,V,lam): return A(q)*(1-1j*V*q)/(1-1j*lam*V*q)

def S0(q,eps):
    # (exp(i q eps)-1 - i q eps)/(i q), regular at 0
    q=np.asarray(q,dtype=complex)
    out=np.empty_like(q)
    small=np.abs(q)<1e-8
    out[~small]=(np.exp(1j*q[~small]*eps)-1-1j*q[~small]*eps)/(1j*q[~small])
    out[small]=1j*q[small]*eps**2/2  # expansion
    return out

def cn(n,eps):
    kn=1j*(n+0.5)
    S0n=(np.exp(1j*kn*eps)-1-1j*kn*eps)/(1j*kn)
    # residue of 1/Mm at kn: (-1)^n/(sqrt(pi) n! i Gamma(1/2-n))
    res=(-1)**n/(np.sqrt(np.pi)*special.factorial(n)*1j*special.gamma(0.5-n))
    return S0n*res

def Gplus(q,eps,N=80):
    # G(q)=S0(q)/Mm(q) + sum c_n/(k_n-q)
    G=S0(q,eps)/Mm(q)
    s=0
    for n in range(N):
        kn=1j*(n+0.5)
        s=s+cn(n,eps)/(kn-q)
    return G+s

def uhat_prime(q,eps,a,V,lam,N=80):
    q=np.asarray(q,dtype=complex)
    # -a/M+ - i q Gplus/M+
    return -a/Mp(q,V,lam) - 1j*q*Gplus(q,eps,N)/Mp(q,V,lam)

# direct Cauchy check of Fplus
def Fplus_direct(q,eps,a,K=300):
    # (1/2pi i) int F(k)/(k-q) dk +? boundary; use Plemelj at q real:
    # F_+ = 1/2 F + (1/(2pi i)) PV int F(k)/(k-q)
    def F(k):
        return ((np.exp(1j*k*eps)-1)+a)/(1j*k*Mm(k))
    val=0.5*F(q)
    # PV integral real axis split
    f=lambda k: F(k)/(k-q)
    # avoid pole by points
    re=integrate.quad(lambda k: f(k).real,-K,q-1e-6,limit=4000,epsabs=1e-10)[0]
    re+=integrate.quad(lambda k: f(k).real,q+1e-6,K,limit=4000,epsabs=1e-10)[0]
    im=integrate.quad(lambda k: f(k).imag,-K,q-1e-6,limit=4000,epsabs=1e-10)[0]
    im+=integrate.quad(lambda k: f(k).imag,q+1e-6,K,limit=4000,epsabs=1e-10)[0]
    val+=(re+1j*im)/(2*np.pi*1j)
    return val

if __name__=='__main__':
    eps=0.125; a=0.4454135
    print('check Fplus at q=0.7,1.5')
    for q in [0.7,1.5]:
        Fp=a/(1j*q)+Gplus(q,eps,N=60)
        Fd=Fplus_direct(q,eps,a,K=200)
        print(q,Fp,Fd,Fp-Fd)
