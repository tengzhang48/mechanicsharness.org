import numpy as np
import eig_g
from eig_g import tquad, tquad as tq
from gbcore import Rfun
# monkey-run eigen_g but check matrices
import scipy.linalg as sla

def build(Nv, Mt):
    def Tbasis(xi, mmax):
        return np.cos(np.outer(np.arccos(np.clip(2*xi-1,-1,1)), np.arange(mmax+1)))
    k = np.arange(1, Nv+2)
    xc = 0.5*(1-np.cos(np.pi*k/(Nv+1)))
    A = np.zeros((Nv+1,Nv+1)); F = np.zeros((Nv+1,Nv+1))
    for i,x in enumerate(xc):
        t, wt = tq(0.0,1.0,Mt); xi = 1-t**2
        Bq = Tbasis(xi,Nv); Bx = Tbasis(np.array([x]),Nv)[0]
        Gq = xi[:,None]/t[:,None]*Bq
        Gx = (x/np.sqrt(1-x))*Bx
        num = t[:,None]*Gq - t[:,None]*Gx[None,:]
        den = x-1+t**2
        I1 = 2*(wt[:,None]*num/den[:,None]).sum(axis=0)
        I2 = Gx*np.log(x/(1-x))
        R = Rfun(x,xi)
        I3 = 2*((wt*R)[:,None]*(t[:,None]*Gq)).sum(axis=0)
        A[i] = I1+I2+I3
        w1 = 2*wt[:,None]*(1-t**2)[:,None]*Bq
        Ig = w1.sum(axis=0); Itg = (w1*xi[:,None]).sum(axis=0)
        tx, wtx = tq(np.sqrt(1-x),1.0,Mt)
        xix = 1-tx**2
        Bqx = Tbasis(xix,Nv)
        Ix = 2*(wtx[:,None]*(1-tx**2)[:,None]*Bqx*((x-xix)[:,None]**2/2)).sum(axis=0)
        F[i] = Itg*x + (-0.5*Ig)*x**2 + Ix
    return A, F, xc

A,F,xc = build(8, 400)
print("A finite:", np.all(np.isfinite(A)), "F finite:", np.all(np.isfinite(F)))
print("A[3]:", A[3][:5])
print("F[3]:", F[3][:5])
