import numpy as np
from numpy.polynomial.chebyshev import chebval
from eig_g import tquad
from gbcore import Rfun, cheb_lobatto_01, bary_weights, lagrange_eval
import scipy.linalg as sla

Nv, Mt = 48, 800
def Tbasis(xi, mmax):
    return np.cos(np.outer(np.arccos(np.clip(2*xi-1,-1,1)), np.arange(mmax+1)))
k = np.arange(1,Nv+2); xc = 0.5*(1-np.cos(np.pi*k/(Nv+2)))
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
lam, V = sla.eig(A,-F)
sel = (np.abs(lam.imag)<1e-8*np.maximum(np.abs(lam.real),1e-3)) & (lam.real>1e-6)
idx = np.where(sel)[0]; order = idx[np.argsort(lam[idx].real)]
Nm = 12

# vectorized f-mode on a dense grid, then barycentric eval
xf = cheb_lobatto_01(240); wf = bary_weights(xf); Lf = lagrange_eval(xf, wf, xf)
t1, wt1 = tquad(0.0,1.0,Mt); xiq = 1-t1**2
def fvals(a):
    vq = chebval(2*xiq-1,a)
    Ig = 2*np.sum(wt1*(1-t1**2)*vq); Itg = 2*np.sum(wt1*(1-t1**2)*xiq*vq)
    out = []
    for x in xf:
        if x == 0.0: out.append(0.0); continue
        tx, wtx = tquad(np.sqrt(1-x),1.0,64); xix = 1-tx**2
        vx = chebval(2*xix-1,a)
        Ix = 2*np.sum(wtx*(1-tx**2)*vx*(x-xix)**2/2)
        out.append(Itg*x - 0.5*Ig*x**2 + Ix)
    return np.array(out)

lams_m = lam[order[:Nm]].real
fgrid = np.array([fvals(np.real(V[:,j])) for j in order[:Nm]])  # (Nm, 241)
def fmode_i(i, xq):
    L = lagrange_eval(xf, wf, xq)
    return L @ fgrid[i]

xg = np.linspace(0.001,0.999,199)
Fm = np.array([fmode_i(n,xg) for n in range(Nm)])
c, *_ = np.linalg.lstsq(Fm.T, np.ones_like(xg), rcond=None)
recon = Fm.T @ c
print("IC reconstruction: ||recon-1||_inf =", np.max(np.abs(recon-1)))
print("c1..c4:", c[:4])
taus = np.linspace(0.3, 1.2, 10)
for x in [0.25,0.5,0.9]:
    vals = np.array([sum(c[n]*fmode_i(n,np.array([x]))[0]*np.exp(-lams_m[n]*t) for n in range(Nm)) for t in taus])
    sl = np.polyfit(taus[2:], np.log(vals[2:]),1)[0]
    print(f"x={x}: late slope = {sl:.4f} (expect {-lams_m[0]:.4f})")
l1 = lams_m[0]; print(f"1/lam1 = {1/l1:.5f}")
def s05(tau): return sum(c[n]*fmode_i(n,np.array([0.5]))[0]*np.exp(-lams_m[n]*tau) for n in range(Nm))
for frac in [0.5,0.05]:
    lo,hi = 0.0,5.0
    for _ in range(60):
        mid=0.5*(lo+hi)
        if s05(mid)>frac: lo=mid
        else: hi=mid
    print(f"s(0.5)={frac}: tau={0.5*(lo+hi):.4f}")
