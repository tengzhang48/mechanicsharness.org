# Independent recomputation of the Hui-Long steady Mode III strip scalars from the stated equations.
import numpy as np
from scipy import integrate, optimize
lam=0.1
def rhs(eps,V,lam=lam):
    # (2/pi^2) int_0^eps [1-(1-lam)e^{-xi/V}] log((1+al)/(1-al)) dxi, xi = eps t^2
    def f(t):
        xi=eps*t*t
        al2=-np.expm1(-(eps-xi))*np.exp(-xi)/(-np.expm1(-eps))   # (e^-xi - e^-eps)/(1-e^-eps)
        al=np.sqrt(max(al2,0.0))
        L=np.log1p(al)-np.log1p(-al)
        return (1-(1-lam)*np.exp(-xi/V))*L*2*eps*t
    return 2/np.pi**2*integrate.quad(f,0,1,limit=400,epsabs=1e-14,epsrel=1e-13)[0]
def eps_of(V,d=0.2,lam=lam):
    return optimize.brentq(lambda e: rhs(e,V,lam)-d, 1e-6, 60, xtol=1e-14)
for V in (0.001,0.1,10):
    e=eps_of(V); a=2/np.pi*np.arccos(np.exp(-e/2)); print('V=%g eps=%.10f a=%.10f R=%.10f'%(V,e,a,a*a/0.2))
d=0.2; print('slow limit eps',-2*np.log(np.cos(np.pi*np.sqrt(d)/2)))
def k(V):  # direct series with tail estimate
    n=np.arange(1,2_000_001); m=2*n-1
    s=np.sum(1.0/(m*(m+2.0/V)))
    return 8/np.pi**2*(s+1/(2*(2*n[-1]+1)))   # tail ~ sum 1/m^2 beyond
def Vstar(d,lam=lam): return optimize.brentq(lambda V: (1-lam)*k(V)-(1-d),0.01,1e4,xtol=1e-12)
for dd in (0.15,0.2,0.3): print('d=%.2f V*=%.10f'%(dd,Vstar(dd)))
# Azure check: relation a vs eps
