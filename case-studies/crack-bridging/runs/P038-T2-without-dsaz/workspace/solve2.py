import numpy as np
from scipy.optimize import brentq

b=0.025; E1=110e9; h1=0.60e-3; E2=190e9; h2=0.90e-3
a0=20.0e-3; sig=0.20e6; dc=0.20e-3; G0=1.40
EI1=E1*b*h1**3/12.0; EI2=E2*b*h2**3/12.0
R=1.0/EI1+1.0/EI2
G=G0+sig*dc/2.0
M=np.sqrt(2.0*b*G/R)
k=b*sig*R/dc; beta=k**0.25

# Model B: front at x=0 has delta=0 (nucleation, intact), delta grows; connector
# first reaches dc at the root x=L  -> boundary to starter cut.
# BCs: x=0: delta=0, delta'=0 ; x=L: delta=dc ; plus moment/shear match with the
# crack/starter-cut region behind: delta''(L)=M R and delta'''(L)=0.
# unknowns A,B,C,D,L with delta=dc+A cosh+b sinh+c cos+d sin.
# x=0: dc + A + C = 0
# x=0: B + D = 0
# x=L: A(cosh-cos)+B(sinh+sin) = 0            (delta=dc)
# x=L: A(sinh+cos?)... ; use delta''(L)=M R and delta'''(L)=0
from scipy.optimize import fsolve
def eqs(v):
    A,B,C,D,L=v
    u=beta*L
    e1=dc+A+C
    e2=B+D
    e3=A*(np.cosh(u)-np.cos(u))+B*(np.sinh(u)+np.sin(u))   # delta(L)-dc=0
    e4=(A*np.cosh(u)+B*np.sinh(u)-C*np.cos(u)-D*np.sin(u))*beta**2 - M*R
    e5=(A*np.sinh(u)+B*np.cosh(u)+C*np.sin(u)-D*np.cos(u))*beta**3
    return [e1,e2,e3,e4,e5]
sol=fsolve(eqs,[1e-4,1e-4,-1e-4,-1e-4,0.012],full_output=False)
A,B,C,D,L=sol
print("Model B: A=%.4e B=%.4e C=%.4e D=%.4e L=%.4f mm  resid=%s"%(A,B,C,D,1e3*L,np.round(eqs(sol),12)))
def d(x):
    ux=beta*x
    return dc+A*np.cosh(ux)+B*np.sinh(ux)+C*np.cos(ux)+D*np.sin(ux)
xs=np.linspace(0,L,200); ds=d(xs)
print("  delta(0)=%.4f um  delta(L/2)=%.4f um  delta(L)=%.4f um  max=%.4f um at %.4f mm"%(
    1e6*d(0),1e6*d(L/2),1e6*d(L),1e6*ds.max(),1e3*xs[np.argmax(ds)]))