import numpy as np
from scipy import special, integrate
def A(q): return np.exp(0.5*np.log(np.pi)+special.loggamma(1-1j*q)-special.loggamma(0.5-1j*q))
def Mm(q): return A(-q)
def I1(K,eps):
    f=lambda k: (np.exp(1j*k*eps)-1)/(1j*k*Mm(k))
    return 2*integrate.quad(lambda k: f(k).real,0,K,limit=60000,epsabs=1e-13,epsrel=1e-12)[0]
def ext(eps,Ks):
    vals=np.array([I1(int(K),eps) for K in Ks])
    M=np.array([[1,-1/np.sqrt(K),-1/K,-K**-1.5] for K in Ks])
    sol=np.linalg.solve(M,vals)
    return sol[0],vals
for eps in [0.05,0.1,0.125,0.15,0.3,0.5,1.0]:
    a3,_=ext(eps,[800,1600,3200])
    a4,v4=ext(eps,[800,1600,3200,6400])
    a5,v5=ext(eps,[1600,3200,6400,12800])
    print('eps',eps,'a3',a3/np.pi,'a4',a4/np.pi,'a5',a5/np.pi,'lastvals',v5)
