import numpy as np, scipy.sparse as spa, scipy.sparse.linalg as spla, fem_j as FJ

nu=0.3;E=1.0;s0=1.0;h=1.0

def run(ell,nx,ny,Xb,Xf):
    K,F,fixed,xs,ys,nxn,nyn,nid,shape,gauss=FJ.build(ell,nx,ny,Xb,Xf)
    ndof=K.shape[0]; free=np.ones(ndof,bool); free[fixed]=False
    u=np.zeros(ndof); u[free]=spla.spsolve(K[free][:,free].tocsc(),F[free])
    return u,xs,ys,nxn,nyn,nid

def J_total(u,xs,ys,nxn,nyn,nid,r,halfheight=None):
    # J = int_Gamma [W n1 - sigma_ij n_j u_i,1] ds, TOTAL stress/sigma,
    # W = total strain energy density (plane stress, foundation included)
    Ds=E/(1-nu**2); k=E*h/(ell_dummy**2*(1-nu**2))
    U1=np.array([[u[2*nid(i,j)] for j in range(nyn)] for i in range(nxn)])
    U2=np.array([[u[2*nid(i,j)+1] for j in range(nyn)] for i in range(nxn)])
    def grad(i,j):
        i1=min(i+1,nxn-1);i0=max(i-1,0);j1=min(j+1,nyn-1);j0=max(j-1,0)
        u1x=(U1[i1,j]-U1[i0,j])/(xs[i1]-xs[i0]);u2x=(U2[i1,j]-U2[i0,j])/(xs[i1]-xs[i0])
        u1y=(U1[i,j1]-U1[i,j0])/(ys[j1]-ys[j0]);u2y=(U2[i,j1]-U2[i,j0])/(ys[j1]-ys[j0])
        e11=u1x;e22=u2y;e12=0.5*(u1y+u2x);tr=e11+e22
        s11=s0+Ds*((1-nu)*e11+nu*tr);s22=s0+Ds*((1-nu)*e22+nu*tr)
        s12=Ds*(1-nu)*e12
        # total strain energy density of the increment field
        W=0.5*Ds*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*tr**2)+0.5*k*(U1[i,j]**2+U2[i,j]**2)
        return (u1x,u2x),(s11,s22,s12),W
    dx=xs[1]-xs[0]; dy=ys[1]-ys[0]
    iR=np.searchsorted(xs,r); iL=np.searchsorted(xs,-r)
    eps=2*dy; jlo=max(1,int(round(eps/dy))); jhi=nyn-1-jlo
    J=0.0
    for j in range(jlo,jhi):   # right edge n=(1,0)
        uu,s,W=grad(iR,j); J+=(W - s[0]*uu[0]-s[2]*uu[1])*dy
    for j in range(jlo,jhi):   # left edge n=(-1,0)
        uu,s,W=grad(iL,j); J+=(-W + s[0]*uu[0]+s[2]*uu[1])*dy
    for i in range(iL,iR):     # top n=(0,1)
        uu,s,W=grad(i,jhi); J+=(-(s[2]*uu[0]+s[1]*uu[1]))*(xs[min(i+1,nxn-1)]-xs[i])
    for i in range(iL,iR):     # bottom n=(0,-1)
        uu,s,W=grad(i,jlo); J+=((s[2]*uu[0]+s[1]*uu[1]))*(xs[min(i+1,nxn-1)]-xs[i])
    return J

ell_dummy=1.0
for nx,ny in [(400,80),(800,160)]:
    u,xs,ys,nxn,nyn,nid=run(1.0,nx,ny,20.0,10.0)
    print(f'nx={nx}')
    for r in [0.25,0.5,1.0,2.0]:
        J=J_total(u,xs,ys,nxn,nyn,nid,r)
        print(f'  r={r} J={J:.6f} J/(s0^2)={J/s0**2:.6f} G*Ebar/(s0^2 ell)={J/(s0**2)}*{1/(1-0.3**2):.4f}={J/(s0**2)/(1-0.3**2):.6f}')