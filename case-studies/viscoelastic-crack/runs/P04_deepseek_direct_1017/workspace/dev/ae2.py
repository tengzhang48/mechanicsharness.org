import numpy as np
from scipy import special, integrate
import time

def A(q): return np.exp(0.5*np.log(np.pi)+special.loggamma(1-1j*q)-special.loggamma(0.5-1j*q))
def Mm(q): return A(-q)

def I2(K=400):
    f=lambda k: 1/(1j*k*Mm(k))+1j/k
    return 2*integrate.quad(lambda k: f(k).real,0,K,limit=4000,epsabs=1e-13,epsrel=1e-13)[0]

def I1(eps,K=400):
    f=lambda k: (np.exp(1j*k*eps)-1)/(1j*k*Mm(k))
    return 2*integrate.quad(lambda k: f(k).real,0,K,limit=4000,epsabs=1e-13,epsrel=1e-13)[0]

if __name__=='__main__':
    t0=time.time(); I2v=I2(); print('I2',I2v,'time',time.time()-t0)
    for eps in [0.02,0.05,0.1,0.136,0.141,0.2,0.3,0.5,1.0,2.0,5.0,10.0]:
        t0=time.time(); v=I1(eps); a=-v/I2v
        print('eps',eps,'I1',v,'a',a,'time',time.time()-t0)
