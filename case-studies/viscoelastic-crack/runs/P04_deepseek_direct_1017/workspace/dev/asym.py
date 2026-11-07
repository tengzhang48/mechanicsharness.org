import numpy as np
from scipy import special
from whfft import Mm,S0
def res1_arr(N):
    r=np.empty(N,dtype=complex); r[0]=-1j/np.pi
    for n in range(1,N): r[n]=r[n-1]*(n-0.5)/n
    return r
def Gplus_naive(q,eps,N):
    r=res1_arr(N); n=np.arange(N); s=n+0.5; kn=1j*s
    S0n=(np.exp(-s*eps)-1+s*eps)/s
    return S0(q,eps)/Mm(q)+np.sum(r*S0n/(kn-q))
def Gplus_resum(q,eps,N=100):
    r=res1_arr(N); n=np.arange(N); s=n+0.5; kn=1j*s
    G=S0(q,eps)/Mm(q)
    return G-eps/Mm(q)+np.sum(r*(np.exp(-s*eps)-1)/(s*(kn-q)))
eps=0.125; a=0.4454135
for q in [1,5,10,50,100,500]:
    naive=Gplus_naive(q,eps,1000000)
    res=Gplus_resum(q,eps,100)
    print('q',q,'naive',naive,'resum',res,'-a/(iq)',-a/(1j*q))
