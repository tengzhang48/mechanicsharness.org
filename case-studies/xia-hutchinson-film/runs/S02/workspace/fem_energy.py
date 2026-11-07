import numpy as np, scipy.sparse as spa, scipy.sparse.linalg as spla
import fem_j as FJ

nu=0.3;E=1.0;s0=1.0;h=1.0;H=1.0

def solve_tip(ell,nx,ny,Xback,Xfront,xtip):
    # crack faces loaded for x1<xtip; far ahead fixed at Xfront
    K,F,fixed,xs,ys,nxn,nyn,nid,shape,gauss=FJ.build(ell,nx,ny,Xback,Xfront)
    # rebuild loads with tip at xtip: FJ.build used tip=0. Redo here:
    F[:]=0.0
    i0=np.searchsorted(xs,xtip)
    for i in range(0,i0):
        Le=xs[i+1]-xs[i]
        F[2*nid(i,0)+1]+= s0*Le/2
        F[2*nid(i+1,0)+1]+= s0*Le/2
        F[2*nid(i,ny)+1]+= -s0*Le/2
        F[2*nid(i+1,ny)+1]+= -s0*Le/2
    ndof=K.shape[0]
    free=np.ones(ndof,bool); free[fixed]=False
    u=np.zeros(ndof)
    u[free]=spla.spsolve(K[free][:,free].tocsc(),F[free])
    # total potential energy Pi = ∫[h s0 tr eps + h/2 Dsig:eps + k/2 u^2]
    Ds=E/(1.0-nu**2)
    k=E*h/(ell**2*(1.0-nu**2))
    U1=np.array([[u[2*nid(i,j)] for j in range(nyn)] for i in range(nxn)])
    U2=np.array([[u[2*nid(i,j)+1] for j in range(nyn)] for i in range(nxn)])
    Pi=0.0
    for i in range(nxn):
        for j in range(nyn):
            i1=min(i+1,nxn-1);i0=max(i-1,0);j1=min(j+1,nyn-1);j0=max(j-1,0)
            u1x=(U1[i1,j]-U1[i0,j])/(xs[i1]-xs[i0]);u2x=(U2[i1,j]-U2[i0,j])/(xs[i1]-xs[i0])
            u1y=(U1[i,j1]-U1[i,j0])/(ys[j1]-ys[j0]);u2y=(U2[i,j1]-U2[i,j0])/(ys[j1]-ys[j0])
            e11=u1x;e22=u2y;e12=0.5*(u1y+u2x);tr=e11+e22
            cell=(xs[i1]-xs[i0])*(ys[j1]-ys[j0]) if False else (xs[1]-xs[0])*(ys[1]-ys[0])
            W=h*s0*tr+0.5*h*Ds*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*tr**2)+0.5*k*(U1[i,j]**2+U2[i,j]**2)
            Pi+=W*cell
    return Pi,xs,ys,nxn,nyn,nid

ell=1.0
nx,ny=400,80;Xback=20.0;Xfront=10.0
Pi0,_,_,_,_,_=solve_tip(ell,nx,ny,Xback,Xfront,0.0)
da=0.05
Pi1,_,_,_,_,_=solve_tip(ell,nx,ny,Xback,Xfront,-da)
# advance tip by da (to -da means crack longer by da) -> release dPi
# per period H, per tip: G*h = -(Pi1-Pi0)/da  (Pi1 has longer crack => lower)
G = -(Pi1-Pi0)/da/h
print(f"Pi0={Pi0:.6f} Pi1={Pi1:.6f}")
print(f"G (per tip, per unit advance) = {G:.6f}")
print(f"G*Ebar/(s0^2 ell) = {G*(1/(1-0.3**2))/(s0**2*ell):.6f}")
print(f"K/(s0 sqrt(ell)) = sqrt(G*E/(s0^2 ell)) = {np.sqrt(G*1.0/(s0**2*ell)):.6f}")