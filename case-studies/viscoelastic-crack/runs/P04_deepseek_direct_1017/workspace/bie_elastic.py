import numpy as np
from numpy.polynomial.legendre import leggauss

def g(u):
    # (1/pi^2) ln|coth(u/4)|, handle u=0
    u=np.asarray(u,dtype=float)
    out=np.zeros_like(u)
    m=np.abs(u)>1e-12
    out[m]=(1/np.pi**2)*np.log(np.abs(1/np.tanh(u[m]/4)))
    return out

def g_log(u):
    return (1/np.pi**2)*np.log(4/np.abs(u))

def int_log(a,b,X):
    # integral_a^b ln(4/|X-s|) ds
    def F(s):
        u=X-s
        return u*np.log(4/np.abs(u))-u
    return F(b)-F(a)

def cell_integral(X,a,b,delta=0.5,nq=8):
    # integral_a^b g(X-s) ds with log subtraction around s=X
    # if interval contains X within delta, split into log part only where |u|<delta
    # simple: use adaptive split at X±delta
    pts=[]
    if a < X-delta < b: pts.append(X-delta)
    if a < X+delta < b: pts.append(X+delta)
    edges=[a]+sorted(pts)+[b]
    tot=0.0
    xg,wg=leggauss(nq)
    for l,r in zip(edges[:-1],edges[1:]):
        if r-l<=0: continue
        mid=(l+r)/2; half=(r-l)/2
        ss=mid+half*xg
        us=X-ss
        val=np.where(np.abs(us)<delta, g_log(us), g(us))
        # where log used, the smooth correction is g-g_log=0 by definition; elsewhere g
        tot += half*np.sum(wg*val)
    return tot

def solve_BIE(eps, L=20.0, N=800, delta=0.5, nq=8):
    edges=np.linspace(0,L,N+1)
    s_mid=0.5*(edges[:-1]+edges[1:])
    X=0.5*(edges[:-1]+edges[1:])  # collocate at cell midpoints
    A=np.zeros((N,N))
    for i,xi in enumerate(X):
        for j in range(N):
            A[i,j]=cell_integral(xi,edges[j],edges[j+1],delta,nq)
    # known negative-side integral: s in [-eps,0]
    b=np.array([-cell_integral(xi,-eps,0,delta,nq) for xi in X])
    # solve least squares
    S,res,rank,sv=np.linalg.lstsq(A,b,rcond=None)
    return edges,s_mid,X,S,sv

if __name__=='__main__':
    for eps in [0.1,0.136,0.2,0.3]:
        edges,s,X,S,sv=solve_BIE(eps,L=15,N=400)
        print('eps',eps,'S(0)',S[0],'S(L)',S[-1],'svmin',sv[-1],'svmax',sv[0])
