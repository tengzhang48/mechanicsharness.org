import numpy as np
import fem as F

nu=0.3; E=1.0; s0=1.0; h=1.0

def fields(u,nxn,nyn,xs,ys,nid):
    # compute sigma increment + total, and gradients at nodes (approx via elements). 
    # We'll just need values along a contour: compute via finite differences of nodal values.
    U1=np.zeros((nxn,nyn)); U2=np.zeros((nxn,nyn))
    for i in range(nxn):
        for j in range(nyn):
            U1[i,j]=u[2*nid(i,j)]; U2[i,j]=u[2*nid(i,j)+1]
    return U1,U2

def J_integral(u,nxn,nyn,xs,ys,nid,ell):
    k=E*h/(ell**2*(1.0-nu**2))
    U1,U2=fields(u,nxn,nyn,xs,ys,nid)
    # element strain energies -> integrate W = h/2 Dsig:eps + k/2 u^2 + h s0 tr eps
    # J = integral over closed contour of (W n1 - sigma_ij n_j u_i,1) ds
    # We'll take a rectangular contour around tip (0,0) inside domain.
    # Use path: far ahead vertical? Simplest: a small rectangle.
    Ds=(E/(1.0-nu**2))
    def eps_sig(i,j):
        # central differences in x1,x2 using grid spacing
        dx=xs[i+1]-xs[i]
        dy=ys[j+1]-ys[j]
        i1=min(i+1,nxn-1); i0=max(i-1,0)
        j1=min(j+1,nyn-1); j0=max(j-1,0)
        u1x=(U1[i1,j]-U1[i0,j])/(xs[i1]-xs[i0])
        u2x=(U2[i1,j]-U2[i0,j])/(xs[i1]-xs[i0])
        u1y=(U1[i,j1]-U1[i,j0])/(ys[j1]-ys[j0])
        u2y=(U2[i,j1]-U2[i,j0])/(ys[j1]-ys[j0])
        e11=u1x; e22=u2y; e12=0.5*(u1y+u2x)
        tr=e11+e22
        s11=s0+Ds*((1-nu)*e11+nu*tr)
        s22=s0+Ds*((1-nu)*e22+nu*tr)
        s12=Ds*(1-nu)*e12
        W = h*s0*tr + 0.5*h*Ds*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*tr**2) + 0.5*k*(U1[i,j]**2+U2[i,j]**2)
        return u1x,u2x,e11,e22,e12,tr,s11,s22,s12,W
    print("J-integral contour computation deferred; need stress on contour.")
    return None

ell=1.0
u,nxn,nyn,xs,ys,nid=F.run(ell)
J_integral(u,nxn,nyn,xs,ys,nid,ell)