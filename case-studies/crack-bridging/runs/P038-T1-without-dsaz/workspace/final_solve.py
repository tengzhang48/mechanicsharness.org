import numpy as np
from scipy.optimize import brentq
import explore
Ep=explore.Ep; a=explore.a; t=explore.t; dc=explore.dc; G0=explore.G0

# monotonic branch: Sigma(L) s.t. G(c=a+L)=G0, using SUPERPOSITION K formula
def Ksuper(c,S): return (S-t)*np.sqrt(np.pi*c)+2*t*np.sqrt(c/np.pi)*np.arcsin(a/c)
def Sig(L):
    c=a+L
    # solve Ksuper(c,S)=sqrt(G0*Ep)
    return brentq(lambda S: Ksuper(c,S)-np.sqrt(G0*Ep), 0.0, 1e11)
def delt(a_x, c, S):   # superposition delta at x
    return (4.0/Ep)*((S-t)*np.sqrt(c*c-a_x*a_x) + t*(a_x - c + np.sqrt(c*c-a_x*a_x)))

Si=np.sqrt(Ep*G0/(np.pi*a))
print("Sigma_i = %.6f MPa" % (Si/1e6))
f=lambda L: delt(a, a+L, Sig(L))-dc
Lstar=brentq(f, 0.0, 1.0)
Sstar=Sig(Lstar)
print("L* = %.6f mm" % (Lstar*1e3))
print("Sigma* = %.6f MPa" % (Sstar/1e6))
print("residual delta(a)=%.10e (target %.10e)" % (delt(a,a+Lstar,Sstar), dc))
print("G=%.10e" % (Ksuper(a+Lstar,Sstar)**2/Ep))

# sensitivity: perturb dc by +-0.3% and see L,Sigma change
for fdc in [0.997,1.003]:
    Lr=brentq(lambda L: delt(a,a+L,Sig(L))-dc*fdc, 0.0,1.0)
    print("dc x%.3f -> L=%.4f mm Sigma=%.4f MPa" % (fdc, Lr*1e3, Sig(Lr)/1e6))