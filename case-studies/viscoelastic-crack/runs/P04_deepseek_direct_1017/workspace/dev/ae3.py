import numpy as np
from scipy import special, integrate
def A(q): return np.exp(0.5*np.log(np.pi)+special.loggamma(1-1j*q)-special.loggamma(0.5-1j*q))
def Mm(q): return A(-q)
def I1(K,eps):
    f=lambda k: (np.exp(1j*k*eps)-1)/(1j*k*Mm(k))
    return 2*integrate.quad(lambda k: f(k).real,0,K,limit=30000,epsabs=1e-12,epsrel=1e-11)[0]
def extrap(eps):
    Ks=[800.,1600.,3200.]
    vals=[I1(int(K),eps) for K in Ks]
    # fit I(K)=Iinf - A K^-1/2 - B K^-1  (tail missing = 4c/sqrtK + ...)
    c=4/np.sqrt(2*np.pi)
    # unknown Iinf, and corrections; use known leading: y=vals - c/sqrtK? careful sign
    # I(K)=Iinf - c/sqrtK - d/K
    # use three points solve
    M=np.array([[1,-1/np.sqrt(K),-1/K] for K in Ks])
    sol=np.linalg.solve(M,np.array(vals))
    return sol[0],vals
if __name__=='__main__':
    for eps in [0.05,0.08,0.1,0.11,0.12,0.125,0.13,0.135,0.14,0.15,0.2,0.3,0.5,1.0]:
        Iinf,vals=extrap(eps)
        print('eps %.3f a %.6f I1 %.6f vals'%(eps,Iinf/np.pi,Iinf),['%.5f'%v for v in vals])
