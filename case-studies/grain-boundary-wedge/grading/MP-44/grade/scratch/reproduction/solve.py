#!/usr/bin/env python3
"""
MP-44  --  stress relaxation by grain-boundary diffusion in a thin film.

Reproduce command:  python3 solve.py

Solves the eigenvalue problem for the normal-traction mode f(Z), Z=z/h,

        lambda f(Z) = R[f](Z) = - int_0^1 S(Z,Xi) f'''(Xi) dXi ,
        S(Z,Xi) = 1/(Z-Xi) - 1/(Z+Xi) - 2 Xi (Z-Xi)/(Z+Xi)^3   (PV at Xi=Z),

obtained from  sigma_t = (E* delta D Omega / 4 pi k T) int_0^h S sigma_zzz d zeta
with tau = t/t0,  t0 = 4 pi k T h^3 /(E* delta D Omega):
        sigma_tau(Z,tau) = int_0^1 S(Z,Xi) sigma_XiXiXi dXi ,
so a mode  sigma = f(Z) exp(-lambda tau)  obeys  lambda f = R[f].

Problem BCs: sigma(0)=0 (free surface); sigma'(h)=sigma''(h)=0 (interface).
They force f(0)=0 and f'(1)=f''(1)=0.  Regularity of f''' at Z=0,1 selects the
weight (Z(1-Z))^{-1/2}, which is also the natural space of the finite Hilbert
transform on (0,1).

Method: weighted-Chebyshev spectral collocation.
  g_k = f'''_k = (Z(1-Z))^{-1/2} T_k(2Z-1), k=0,1,2,...
  * Cauchy part closed form:
        PV int_0^1 (Xi(1-Xi))^{-1/2} T_k(2Xi-1)/(Z-Xi) dXi
            = -2 pi U_{k-1}(2Z-1)  (k>=1),  = 0 (k=0)
  * regular part Rreg=(Xi^2-4 Z Xi-Z^2)/(Z+Xi)^3 by paneled Gauss quadrature
    (split at Xi=Z; endpoint sqrt removed by Xi=Z sin^2 phi / Xi=1-(1-Z)sin^2 phi).
  * f_k = triple integral of g_k with f(0)=f'(1)=f''(1)=0:
        f_k(Z)=int_0^Z (Z-s)^2/2 g_k ds + (I0-I1) Z - (I0/2) Z^2,
        I0=int_0^1 g_k ds, I1=int_0^1 (1-s) g_k ds   (s=sin^2 th).
  Collocation at interior Chebyshev points Z_i: R[i,j]=R[f_j](Z_i), M[i,j]=f_j(Z_i).
"""
import numpy as np, json
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1.0,1.0)))
def chebU(k,x):
    x=np.clip(x,-1.0,1.0); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)

xth,wth   = leggauss(3000)     # theta grid (g ds = 2 T_k dtheta)
xphi,wphi = leggauss(300)      # phi grid for the regular kernel part

def Rf(k,Z):
    """R[f_k](Z) = -[ PV int_0^1 g_k/(Z-Xi) dXi + int_0^1 Rreg g_k dXi ]."""
    x=2.0*Z-1.0
    SP = 0.0 if k==0 else -2.0*np.pi*chebU(k-1,x)      # exact Cauchy part
    ph=0.25*np.pi*(xphi+1.0); w=0.25*np.pi*wphi; s=np.sin(ph); c=np.cos(ph)
    X1=Z*s*s; dX1=2.0*Z*s*c                            # Xi in (0,Z)
    g1=chebT(k,2.0*X1-1.0)/np.sqrt(X1*(1.0-X1))
    RP1=np.sum(w*((X1**2-4.0*Z*X1-Z**2)/(Z+X1)**3)*g1*dX1)
    X2=1.0-(1.0-Z)*s*s; dX2=2.0*(1.0-Z)*s*c           # Xi in (Z,1)
    g2=chebT(k,2.0*X2-1.0)/np.sqrt(X2*(1.0-X2))
    RP2=np.sum(w*((X2**2-4.0*Z*X2-Z**2)/(Z+X2)**3)*g2*dX2)
    return -(SP+RP1+RP2)

def theta_int(k,Z,hfun):
    """int_0^Z hfun(s) g_k(s) ds  with  s=sin^2(theta); g_k ds = 2 T_k dtheta."""
    thZ=np.arcsin(np.sqrt(Z)); th=0.5*thZ*(xth+1.0); w=0.5*thZ*wth; s=np.sin(th)**2
    return 2.0*np.sum(w*hfun(s)*chebT(k,np.cos(2.0*th)))

def fk(k,Z):
    """f_k(Z)=int_0^Z (Z-s)^2/2 g_k ds + (I0-I1) Z - (I0/2) Z^2."""
    F2=theta_int(k,Z,lambda s:0.5*(Z-s)**2)
    I0=theta_int(k,1.0,lambda s:np.ones_like(s))
    I1=theta_int(k,1.0,lambda s:(1.0-s))
    return F2+(I0-I1)*Z+(-I0/2.0)*Z**2

def solve(N):
    M=N; th=(np.arange(1,M+1))*np.pi/(M+1); Zc=0.5*(1-np.cos(th))
    Rm=np.zeros((M,M)); Mm=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(N):
            Rm[i,j]=Rf(j,Z); Mm[i,j]=fk(j,Z)
    lam,vec=eig(Rm,Mm)
    return lam,vec,Zc
# --------------------------------------------------------------------- checks
def _chebT_arr(k,y):
    T0=np.ones_like(y); T1=y.copy()
    if k==0: return T0
    if k==1: return T1
    for n in range(2,k+1): T0,T1=T1,2.0*y*T1-T0
    return T1
_xb,_wb=leggauss(800); _xi=0.5*(_xb+1.0); _wt=0.5*_wb; _y=2.0*_xi-1.0
def Rf_brute(k,Z):
    """Independent R[f_k](Z): PV by singularity subtraction, no analytic identity."""
    gv=_chebT_arr(k,_y)/np.sqrt(_xi*(1.0-_xi))
    gZ=_chebT_arr(k,np.array([2.0*Z-1.0]))[0]/np.sqrt(Z*(1.0-Z))
    pv=np.sum(_wt*(gv-gZ)/(Z-_xi))+gZ*np.log(abs(Z/(1.0-Z)))
    reg=np.sum(_wt*((_xi**2-4.0*Z*_xi-Z**2)/(Z+_xi)**3)*gv)
    return -(pv+reg)

def Rf_indep(k,Z):
    """Independent R[f_k](Z): Cauchy part with numerical singularity subtraction
    (NO analytic Hilbert identity), paneled theta-free quadrature at Xi=Z."""
    ph=0.25*np.pi*(xphi+1.0); w=0.25*np.pi*wphi; s=np.sin(ph); c=np.cos(ph)
    gZ=chebT(k,2.0*Z-1.0)/np.sqrt(Z*(1.0-Z))
    X1=Z*s*s; dX1=2.0*Z*s*c
    X2=1.0-(1.0-Z)*s*s; dX2=2.0*(1.0-Z)*s*c
    g1=chebT(k,2.0*X1-1.0)/np.sqrt(X1*(1.0-X1))
    g2=chebT(k,2.0*X2-1.0)/np.sqrt(X2*(1.0-X2))
    pv=(np.sum(w*(g1-gZ)/(Z-X1)*dX1)+np.sum(w*(g2-gZ)/(Z-X2)*dX2)
        +gZ*np.log(abs(Z/(1.0-Z))))
    RP1=np.sum(w*((X1**2-4.0*Z*X1-Z**2)/(Z+X1)**3)*g1*dX1)
    RP2=np.sum(w*((X2**2-4.0*Z*X2-Z**2)/(Z+X2)**3)*g2*dX2)
    return -(pv+RP1+RP2)

if __name__=="__main__":
    N=64
    lam,vec,Zc=solve(N)
    cand=sorted([(l.real,i) for i,l in enumerate(lam)
                 if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    lam1,lam2=cand[0][0],cand[1][0]
    c=vec[:,cand[0][1]].real
    f1=sum(c[j]*fk(j,1.0) for j in range(N))
    half=sum(c[j]*fk(j,0.5) for j in range(N))/f1
    Kinf=1.1215                 # edge crack of length h : K=1.1215 sigma0 sqrt(pi h)
    relax=1.0/lam1

    relmax=0.0; relind=0.0
    for k in [0,1,2,3,5]:
        for Z in [0.15,0.2,0.35,0.5,0.7,0.85]:
            a=Rf(k,Z); b=Rf_brute(k,Z); e=Rf_indep(k,Z)
            relmax=max(relmax,abs(a-b)/max(1e-30,abs(b)))
            relind=max(relind,abs(a-e)/max(1e-30,abs(e)))
    Zt=np.linspace(0.05,0.95,19)
    fv=np.array([sum(c[j]*fk(j,Z) for j in range(N)) for Z in Zt])
    Rv=np.array([sum(c[j]*Rf(j,Z) for j in range(N)) for Z in Zt])
    resid=np.max(np.abs(Rv/fv-lam1))
    cn=np.abs(c)/np.max(np.abs(c))
    shp=np.array([sum(c[j]*fk(j,Z) for j in range(N))/f1 for Z in Zt])

    print("VERIFY R[f] vs 800-pt direct PV integration      : max rel err = %.2e"%relmax)
    print("VERIFY R[f] vs independent paneled PV (no identity): max rel err = %.2e"%relind)
    print("VERIFY eigenmode residual max|R[f]/f - lambda1|              = %.2e"%resid)
    print("VERIFY mode coefficients |c_k|/max (first 12):", np.round(cn[:12],3))
    print("VERIFY f/f(1) at Z:", np.round(Zt,3))
    print("            f/f(1)   :", np.round(shp,4))
    print()
    out={"lambda1":round(float(lam1),4),
         "lambda2":round(float(lam2),3),
         "mode1_shape_half":round(float(half),4),
         "K_inf_norm":Kinf,
         "relaxation_time_t0":round(float(relax),5),
         "method":("weighted-Chebyshev spectral collocation for the nonlocal traction "
                   "eigenproblem lambda f = -int_0^1 S f''' dXi (f''' in "
                   "(Z(1-Z))^{-1/2}T_k(2Z-1), PV Cauchy part analytic, regular part by "
                   "paneled Gauss quadrature)")}
    print(json.dumps(out,indent=2))