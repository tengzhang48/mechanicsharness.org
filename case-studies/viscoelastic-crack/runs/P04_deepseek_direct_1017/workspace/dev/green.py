import numpy as np
from scipy import integrate

def g_el(X):
    X=np.asarray(X,float); out=np.empty_like(X)
    ax=np.abs(X)
    # log|coth(x/4)|, stable
    out=np.where(ax<1e-8, np.log(4.0/np.maximum(ax,1e-300)),
                 np.log(1.0/np.tanh(ax/4.0)))
    return out/np.pi**2

def g_ve_scalar(X,V,lam):
    # lambda*g_el(X)+(1-lam)*(1/V) int_0^inf g_el(X+s) e^{-s/V} ds
    # split to handle log singular at s=-X for X<0
    f=lambda s: g_el(np.array([X+s]))[0]*np.exp(-s/V)/V
    if X>=0:
        val=integrate.quad(f,0,np.inf,epsabs=1e-13,epsrel=1e-13,limit=400)[0]
    else:
        # singular at s=-X>0
        s0=-X
        val=integrate.quad(f,0,s0,epsabs=1e-13,epsrel=1e-13,limit=400,points=[s0])[0]
        val+=integrate.quad(f,s0,np.inf,epsabs=1e-13,epsrel=1e-13,limit=400)[0]
    return lam*g_el(np.array([X]))[0]+(1-lam)*val

if __name__=='__main__':
    lam=0.1
    for V in [0.001,0.1,10.0]:
        print('V',V)
        for X in [0.0,0.5,1.0,2.0,5.0,-0.5,-1.0,-2.0]:
            print('  X',X,'g_ve',g_ve_scalar(X,V,lam),'g_el',g_el(np.array([X]))[0])
