import numpy as np
from scipy.optimize import brentq
import solve as S

Ep=S.Ep; a=S.a; t=S.t; K0=S.K0; dc=S.dc

def solve_SIE(Sr, c, M):
    xa=a/c
    k=np.arange(1,M+1); xi=np.cos(k*np.pi/(M+1)); s2=np.sin(k*np.pi/(M+1))**2
    g=np.where(np.abs(xi)>xa, Sr-t, Sr); F=-(4.0*c/Ep)*g
    b=np.zeros(M+1); U0=np.ones(M); U1=2*xi
    b[1]=(2.0/(M+1))*np.dot(s2*F,U0); b[2]=(2.0/(M+1))*np.dot(s2*F,U1)
    for n in range(3,M+1):
        U2=2*xi*U1-U0; b[n]=(2.0/(M+1))*np.dot(s2*F,U2); U0,U1=U1,U2
    K=-(Ep/4)*np.sqrt(np.pi/c)*b[1:].sum()
    def delt(xi0):
        s=0.0; V0=1.0; V1=2*xi0; s=b[1]*V0
        for n in range(2,M+1):
            s+=b[n]*V1/n
            V2=2*xi0*V1-V0; V0,V1=V1,V2
        return -np.sqrt(1-xi0*xi0)*s
    return K, delt

c_s=S.c_s; S_s=S.S_s
# convergence of spectral at the closed-form solution
for M in [1024, 4096, 8192]:
    K,delt = solve_SIE(S_s, c_s, M)
    print(f"M={M:5d}: K_tip={K/1e6:.5f} (target {K0/1e6:.5f}), delta(a)={delt(a/c_s)*1e6:.5f} um (target 9.0)")

# fully independent 2x2 solve with spectral method only
M=4096
def Sig_from_K(c):
    Kan = lambda Sr: solve_SIE(Sr,c,M)[0]
    return brentq(lambda Sr: Kan(Sr)-K0, 1e6, 159e6, xtol=1e3)
def resid(c):
    Sr = Sig_from_K(c)
    return solve_SIE(Sr,c,M)[1](a/c) - dc
c_spec = brentq(resid, 2.5e-3, 9e-3, xtol=1e-10, rtol=1e-12)
S_spec = Sig_from_K(c_spec)
print("\nspectral-only 2x2 solve: c =", c_spec*1e3, "mm, L =", (c_spec-a)*1e3, "mm, Sigma =", S_spec/1e6, "MPa")
print("closed-form solution   : c =", c_s*1e3, "mm, L =", (c_s-a)*1e3, "mm, Sigma =", S_s/1e6, "MPa")
print("relative differences: dL/L =", abs(c_spec-c_s)/(c_s-a), ", dSigma/Sigma =", abs(S_spec-S_s)/S_s)

# monotone opening profile over the whole crack (normalized)
K,delt = solve_SIE(S_s, c_s, 4096)
xis = np.linspace(0.0, 1.0, 21)
prof = [delt(x*0.99999 if x>=1 else x) for x in xis]
print("\nprofile delta(x)/delta(a), xi=0(center)..1(tip):")
print(np.round(np.array(prof)/delt(a/c_s),4))
