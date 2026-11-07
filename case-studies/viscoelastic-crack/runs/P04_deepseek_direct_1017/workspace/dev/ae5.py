import numpy as np
from scipy import special, integrate
def A(q): return np.exp(0.5*np.log(np.pi)+special.loggamma(1-1j*q)-special.loggamma(0.5-1j*q))
def Mm(q): return A(-q)
def I1(K,eps):
    f=lambda k: (np.exp(1j*k*eps)-1)/(1j*k*Mm(k))
    return 2*integrate.quad(lambda k: f(k).real,0,K,limit=100000,epsabs=1e-13,epsrel=1e-12)[0]
def ext(eps,Ks,pows):
    vals=np.array([I1(int(K),eps) for K in Ks])
    M=np.array([[K**(-p) for p in pows] for K in Ks])
    sol=np.linalg.solve(M,vals)
    return sol[0],vals
for eps in [0.05,0.1,0.125,0.15,0.3,0.5,1.0]:
    a3,_=ext(eps,[800,1600,3200],[0,0.5,1.0])
    a4,v4=ext(eps,[800,1600,3200,6400],[0,0.5,1.0,1.5])
    a5,v5=ext(eps,[1600,3200,6400,12800],[0,0.5,1.0,1.5])
    a6,v6=ext(eps,[1600,3200,6400,12800],[0,0.5,1.0,1.5,2.0][:4]) if False else (0,0)
    print('eps',eps,'a3 %.7f a4 %.7f a5 %.7f'%(a3/np.pi,a4/np.pi,a5/np.pi))
    print('   vals',v5)
