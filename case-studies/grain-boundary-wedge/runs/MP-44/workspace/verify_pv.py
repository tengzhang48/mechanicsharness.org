import numpy as np
from numpy.polynomial.legendre import leggauss

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1,1)))
def chebU(k,x):
    x=np.clip(x,-1,1)
    if k==0: return np.ones_like(x)
    th=np.arccos(x); s=np.sin(th)
    out=np.where(s>1e-12, np.sin((k+1)*th)/np.where(s>1e-12,s,1.0), chebU_limit(k,x))
    return out
def chebU_limit(k,x):
    # limit at x=+-1: U_k(1)=k+1, U_k(-1)=(-1)^k (k+1)
    return np.where(x>0, k+1.0, ((-1)**k)*(k+1.0))*np.ones_like(x)

# verify PV identities by direct numerical PV with singularity subtraction
xg,wg=leggauss(4000); xgp=xg  # x in (-1,1)
def pv_direct(k,fam,y):
    # compute PV int_{-1}^1 base_k(y')/(y-y') dy' where base includes weight
    if fam=='A':
        base=lambda t:(1-t**2)**(-0.5)*chebT(k,t)
    else:
        base=lambda t:(1-t**2)**0.5*chebU(k,t)
    # PV int base(t)/(y-t) dt = int (base(t)-base(y))/(y-t)dt + base(y) PV int dt/(y-t)
    num=(base(xgp)-base(y)); den=(y-xgp)
    reg=np.sum(wg*num/den)
    pvlog=base(y)*np.log(np.abs((y+1)/(y-1)))  # PV int_{-1}^1 dt/(y-t)= ln|(y+1)/(y-1)|
    return reg+pvlog

for k in [1,2,3]:
    for y in [-0.5,0.0,0.5]:
        dA=pv_direct(k,'A',y); idA=-np.pi*chebU(k-1,y)
        print("A k=%d y=%.1f direct=%.6f ident=%.6f"%(k,y,dA,idA))
for k in [0,1,2]:
    for y in [-0.5,0.0,0.5]:
        dB=pv_direct(k,'B',y); idB=np.pi*chebT(k+1,y)
        print("B k=%d y=%.1f direct=%.6f ident=%.6f"%(k,y,dB,idB))