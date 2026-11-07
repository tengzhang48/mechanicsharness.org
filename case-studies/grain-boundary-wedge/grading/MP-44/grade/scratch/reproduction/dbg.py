import numpy as np
from numpy.polynomial.legendre import leggauss

def chebT(k,x): return np.cos(k*np.arccos(np.clip(x,-1,1)))
def chebU(k,x):
    x=np.clip(x,-1,1); th=np.arccos(x); s=np.sin(th)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(s>1e-10, np.sin((k+1)*th)/np.where(s>1e-10,s,1.0), 0.0)
def gB(k,Z):
    x=2*Z-1; return (1-x**2)**0.5*chebU(k,x)
def gA(k,Z):
    x=2*Z-1; return (1-x**2)**(-0.5)*chebT(k,x)

# x,y quadratures
xgx,wgx=leggauss(400)
# Xi quadrature for arbitrary smooth-ish g (use substitution Xi=sin^2 to tame sqrt)
def brute_parts(gfun,Z):
    # tanh-sinh-ish: use Gauss on (0,1) with many pts
    xg2,wg2=leggauss(4000); xg2=0.5*(xg2+1); wg2=0.5*wg2
    gv=gfun(xg2)
    d=np.where(np.abs(Z-xg2)>1e-14, Z-xg2, 1.0)
    sing=np.sum(wg2*(gv-gfun(Z))/d)+gfun(Z)*np.log(Z/(1-Z))
    reg=np.sum(wg2*(xg2**2-4*Z*xg2-Z**2)/(Z+xg2)**3*gv)
    return sing,reg

def id_parts(k,fam,Z):
    x=2*Z-1
    if fam=='A': SP=0.0 if k==0 else np.pi*chebU(k-1,x)
    else: SP=-np.pi*chebT(k+1,x)
    y=xgx; g=gA(k,(y+1)/2) if fam=='A' else gB(k,(y+1)/2)
    RP=np.sum(wgx*(-1.0/(x+y+2)-2*(y+1)*(x-y)/(x+y+2)**3)*g)
    return SP,RP

for fam,k in [('B',2),('A',1)]:
    for Z in [0.2,0.5,0.7]:
        s1,r1=brute_parts(lambda W:gB(k,W) if fam=='B' else gA(k,W),Z)
        s2,r2=id_parts(k,fam,Z)
        print(fam,k,round(Z,2),"sing: brute=%.5f ident=%.5f | reg: brute=%.5f ident=%.5f"%(s1,s2,r1,r2))