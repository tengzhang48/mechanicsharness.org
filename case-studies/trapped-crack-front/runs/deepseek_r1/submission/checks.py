import numpy as np
from scipy.optimize import brentq

nu=0.3; K2r=0.6; K3r=0.8; q=0.5
Gcoef=(1-nu**2)*(1+K2r**2)+(1+nu)*K3r**2
w1=(1-nu**2)/Gcoef; w2=(1-nu**2)*K2r**2/Gcoef; w3=(1+nu)*K3r**2/Gcoef
beta=2*(-(4/3)*w2+0.75*w3)

# Part A direct check
C=np.array([-np.pi, -K2r*np.pi, -K3r*np.pi])
S=np.array([0.0, 2*np.pi*K3r, -2*np.pi*K2r])
# relative perturbations = (A/lam)*C_alpha/(K_alpha/K1)
ratio=np.array([1.0,K2r,K3r])
ss=np.linspace(0,1,8,endpoint=False)
eps=1e-6
for s in ss[:4]:
    th=2*np.pi*s
    d=(eps*np.array([C[0]*np.cos(th)+S[0]*np.sin(th),
                     C[1]*np.cos(th)+S[1]*np.sin(th),
                     C[2]*np.cos(th)+S[2]*np.sin(th)])/ratio)
    dG=2*(w1*d[0]+w2*d[1]+w3*d[2])
    pred=eps*(-2*np.pi*np.cos(th)-2*np.pi*beta*np.sin(th))
    print("s=%.3f dG=%.6e pred=%.6e diff=%.2e"%(s,dG,pred,dG-pred))

# Part B check G=Gamma for a finite eps using the first-order kappas
def D(s):
    t=2*np.pi*s; return 1-2*q*np.cos(t)+q*q
def phi(s):
    t=2*np.pi*s; return np.arctan2(-q*np.sin(t),1-q*np.cos(t))
def g(s):
    t=2*np.pi*s; return 2*q*(np.cos(t)-q)/D(s)
def f(s):
    return (0.5*np.log(D(s))-beta*phi(s))/(np.pi*(1+beta**2))
def fp(s):
    t=2*np.pi*s
    return (2*q/(1+beta**2))*(np.sin(t)-beta*(q-np.cos(t)))/D(s)
def kap(s):
    k1=(g(s)-beta*fp(s))/2
    return k1, K2r*k1-K3r*fp(s), K3r*k1+K2r*fp(s)
for eps in [1e-2,1e-3,1e-4]:
    worst=0
    for s in np.linspace(0,1,101)[:-1]:
        k1,k2,k3=kap(s)
        d=np.array([eps*k1,eps*k2/K2r,eps*k3/K3r])
        dG=2*(w1*d[0]+w2*d[1]+w3*d[2])
        worst=max(worst,abs(dG-eps*g(s)))
    print("eps",eps,"max residual",worst)

# Verify the front shape solves the equilibrium equation using the actual
# first-order relation kappa1 = -H[f']/2, with f' and H[f'] from FFT.
M=4096
s=np.arange(M)/M
fv=f(s); fpr=fp(s); gv=g(s)
F=np.fft.rfft(fv); n=np.arange(len(F))
Hf=np.fft.irfft(2*np.pi*n*F,n=M)
print("max|-Hf/2+beta fp/2-g/2| =",np.max(np.abs(-Hf/2+beta*fpr/2-gv/2)))

# kappa definitions from equilibrium vs direct -H f'/2
k1dir=-Hf/2
k1eq=(gv-beta*fpr)/2
print("max|k1dir-k1eq| =",np.max(np.abs(k1dir-k1eq)))

# Check Fourier coefficients of f against closed form
fhat_closed=np.zeros(len(F),dtype=complex)
for k in range(1,len(F)):
    fhat_closed[k]=M*q**k*(-1-1j*beta)/(2*np.pi*k*(1+beta**2))
f_recon=np.fft.irfft(fhat_closed,n=M)
print("max f Fourier recon error =",np.max(np.abs(f_recon-fv)))
