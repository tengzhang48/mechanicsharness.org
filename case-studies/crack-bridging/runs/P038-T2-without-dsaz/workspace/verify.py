import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

b=0.025; E1=110e9; h1=0.60e-3; E2=190e9; h2=0.90e-3
sig=0.20e6; dc=0.20e-3; G0=1.40
EI1=E1*b*h1**3/12.0; EI2=E2*b*h2**3/12.0
R=1.0/EI1+1.0/EI2
tot=G0+sig*dc/2.0
M=np.sqrt(2.0*b*tot/R)
k=b*sig*R/dc; beta=k**0.25
A=M*R/(2*beta**2); C=-A
def resid(l):
    u=beta*l;B=-A*(np.sinh(u)+np.sin(u))/(np.cosh(u)+np.cos(u))
    return A*(np.cosh(u)-np.cos(u))+B*(np.sinh(u)+np.sin(u))+dc
ls=np.linspace(1e-4,0.1,8000);v=[resid(l) for l in ls];sc=np.where(np.diff(np.sign(v)))[0]
ell=brentq(resid,ls[sc[0]],ls[sc[0]+1]);u=beta*ell
B=-A*(np.sinh(u)+np.sin(u))/(np.cosh(u)+np.cos(u));D=B
def d(x):
    ux=beta*x; return dc+A*np.cosh(ux)+B*np.sinh(ux)+C*np.cos(ux)+D*np.sin(ux)
def d1(x):
    ux=beta*x; return beta*(A*np.sinh(ux)+B*np.cosh(ux)-C*np.sin(ux)+D*np.cos(ux))
def d2(x):
    ux=beta*x; return beta**2*(A*np.cosh(ux)+B*np.sinh(ux)-C*np.cos(ux)-D*np.sin(ux))
def d3(x):
    ux=beta*x; return beta**3*(A*np.sinh(ux)+B*np.cosh(ux)+C*np.sin(ux)-D*np.cos(ux))
# BCs
print('BC tip  delta(0)=dc? ', 1e6*d(0), 'dc_um=',1e6*dc)
print('BC tip  delta\'\'(0)=M R? ', d2(0), M*R)
print('BC tip  delta\'\'\'(0)=0? ', d3(0))
print('BC root delta(L)=0? ', 1e6*d(ell))
print('BC root delta\'(L)=0? ', d1(ell))
# energy: cohesive dissipation along the zone (per unit width) and area
work_per_area = sig*dc/2
print('cohesive work per area = %.4f J/m2'%work_per_area)
# integral of t*|d_delta/dx|*dx over zone = area under t-delta = 1/2 sig dc (for each point)
# total dissipated per unit width = integral t(delta) * (d(delta)/dx)?? no; it's per unit area = 0.5 sig dc
print('total G supply = %.4f J/m2'%tot)
print('ANSWER: M=%.6f N.m  L=%.4f mm  d(L/2)=%.4f um'%(M,1e3*ell,1e6*d(ell/2)))
