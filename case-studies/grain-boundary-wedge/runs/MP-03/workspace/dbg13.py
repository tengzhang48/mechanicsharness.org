import numpy as np
from numpy.polynomial.chebyshev import chebval
from eig_g import tquad
from gbcore import Rfun
import scipy.linalg as sla

def build_AFP(Nv, Mt, xc):
    def Tbasis(xi, mmax):
        return np.cos(np.outer(np.arccos(np.clip(2*xi-1,-1,1)), np.arange(mmax+1)))
    A = np.zeros((Nv+1,Nv+1)); F = np.zeros((Nv+1,Nv+1))
    for i,x in enumerate(xc):
        t, wt = tquad(0.0,1.0,Mt); xi = 1-t**2
        Bq = Tbasis(xi,Nv); Bx = Tbasis(np.array([x]),Nv)[0]
        Gq = xi[:,None]/t[:,None]*Bq; Gx = (x/np.sqrt(1-x))*Bx
        den = x-1+t**2
        I1 = 2*(wt[:,None]*(t[:,None]*Gq - t[:,None]*Gx[None,:])/den[:,None]).sum(axis=0)
        I2 = Gx*np.log(x/(1-x))
        I3 = 2*((wt*Rfun(x,xi))[:,None]*(t[:,None]*Gq)).sum(axis=0)
        A[i] = I1+I2+I3
        w1 = 2*wt[:,None]*(1-t**2)[:,None]*Bq
        Ig = w1.sum(axis=0); Itg = (w1*xi[:,None]).sum(axis=0)
        tx, wtx = tquad(np.sqrt(1-x),1.0,Mt); xix = 1-tx**2; Bqx = Tbasis(xix,Nv)
        Ix = 2*(wtx[:,None]*(1-tx**2)[:,None]*Bqx*((x-xix)[:,None]**2/2)).sum(axis=0)
        F[i] = Itg*x + (-0.5*Ig)*x**2 + Ix
    return A,F

def run(Nv, Mt, xctype):
    if xctype=='cheb':
        k = np.arange(1,Nv+2); xc = 0.5*(1-np.cos(np.pi*k/(Nv+2)))
    else:
        xc = np.linspace(1.0/(Nv+1), 1.0-1.0/(Nv+1), Nv+1)
    A,F = build_AFP(Nv,Mt,xc)
    lam,V = sla.eig(A,-F)
    sel = (np.abs(lam.imag)<1e-8*np.maximum(np.abs(lam.real),1e-3)) & (lam.real>1e-6)
    return np.sort(lam[sel].real)

for tag, Nv, Mt, xct in [("cheb Nv=16 Mt=200",16,200,'cheb'),
                         ("cheb Nv=48 Mt=200",48,200,'cheb'),
                         ("cheb Nv=48 Mt=1600",48,1600,'cheb'),
                         ("unif Nv=48 Mt=800",48,800,'unif'),
                         ("unif Nv=32 Mt=800",32,800,'unif')]:
    l = run(Nv,Mt,xct)
    print(f"{tag}: lam1={l[0]:.7f} lam2={l[1]:.5f} lam3={l[2]:.4f}")
