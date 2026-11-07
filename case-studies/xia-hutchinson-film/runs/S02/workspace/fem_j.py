import numpy as np
import scipy.sparse as spa
import scipy.sparse.linalg as spla

# ======================================================================
# CASE A, FEM with J-integral around the tip. Ground truth.
# ======================================================================
nu=0.3; E=1.0; s0=1.0; h=1.0; H=1.0

def build(ell, nx, ny, Xback, Xfront, grad=True):
    k=E*h/(ell**2*(1.0-nu**2))
    xs=np.linspace(-Xback,Xfront,nx+1); ys=np.linspace(0.0,H,ny+1)
    nxn=len(xs); nyn=len(ys)
    nid=lambda i,j:i*nyn+j
    nn=nxn*nyn; ndof=2*nn
    D=(E/(1.0-nu**2))*np.array([[1.0,nu,0.0],[nu,1.0,0.0],[0.0,0.0,(1.0-nu)/2.0]])
    gauss=[(-1/np.sqrt(3),-1/np.sqrt(3)),(1/np.sqrt(3),-1/np.sqrt(3)),
           (1/np.sqrt(3),1/np.sqrt(3)),(-1/np.sqrt(3),1/np.sqrt(3))]
    def shape(xi,eta):
        N=np.array([(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4])
        dN=np.array([[-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4],
                     [-(1-xi)/4,-(1+xi)/4,(1+xi)/4,(1-xi)/4]])
        return N,dN
    rows=[];cols=[];vals=[];F=np.zeros(ndof)
    for i in range(nx):
        for j in range(ny):
            nodes=[nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)]
            xy=np.array([[xs[i],ys[j]],[xs[i+1],ys[j]],[xs[i+1],ys[j+1]],[xs[i],ys[j+1]]])
            ke=np.zeros((8,8))
            for xi,eta in gauss:
                N,dN=shape(xi,eta)
                J=dN@xy; detJ=np.linalg.det(J)
                dNxy=np.linalg.solve(J,dN)
                B=np.zeros((3,8))
                B[0,0::2]=dNxy[0]; B[1,1::2]=dNxy[1]
                B[2,0::2]=dNxy[1]; B[2,1::2]=dNxy[0]
                ke+=(B.T@D@B)*detJ*h
                Nmat=np.zeros((2,8)); Nmat[0,0::2]=N; Nmat[1,1::2]=N
                ke+=(Nmat.T@Nmat)*detJ*k
            for a_ in range(8):
                for b_ in range(8):
                    rows.append(2*nodes[a_//2]+(a_%2));cols.append(2*nodes[b_//2]+(b_%2));vals.append(ke[a_,b_])
    K=spa.csr_matrix((vals,(rows,cols)),shape=(ndof,ndof))
    fixed=[]
    for j in range(nyn):
        n=nid(nx-1,j); fixed+=[2*n,2*n+1]
    # also fix far behind to zero? No: foundation handles it. But to avoid
    # rigid modes, u is pinned by foundation. Keep only far-ahead fixed.
    i0=np.searchsorted(xs,0.0)
    for i in range(0,i0):
        Le=xs[i+1]-xs[i]
        F[2*nid(i,0)+1]+= s0*Le/2
        F[2*nid(i+1,0)+1]+= s0*Le/2
        F[2*nid(i,ny)+1]+= -s0*Le/2
        F[2*nid(i+1,ny)+1]+= -s0*Le/2
    return K,F,fixed,xs,ys,nxn,nyn,nid,shape,gauss

def solve(ell,nx,ny,Xback,Xfront):
    K,F,fixed,xs,ys,nxn,nyn,nid,shape,gauss=build(ell,nx,ny,Xback,Xfront)
    ndof=K.shape[0]
    free=np.ones(ndof,bool); free[fixed]=False
    u=np.zeros(ndof)
    u[free]=spla.spsolve(K[free][:,free].tocsc(),F[free])
    return u,xs,ys,nxn,nyn,nid,shape,gauss

def Jint(ell,u,xs,ys,nxn,nyn,nid,shape,gauss,r=0.5):
    # rectangular contour around tip (0,0): path in x1 in [-r? no,
    # crack occupies x1<0 at x2=0 and x2=H]. A contour from x2=0 to H
    # ahead of tip, back behind at x2=H ... standard J path:
    # choose rectangle: x1=+r (right edge), x1=-r (left edge), and follow
    # around. Simpler: J = integral along a path from bottom-far to top-far
    # enclosing the tip. Use: vertical segment at x1=+r from x2=0..H,
    # horizontal at x2=0 from x1=+r to -r, then x1=-r up, then x2=H back.
    # For a semi-infinite crack under remote load, J is path independent.
    k=E*h/(ell**2*(1.0-nu**2))
    def nodal_vals(Uarr):
        pass
    U1=np.array([[u[2*nid(i,j)] for j in range(nyn)] for i in range(nxn)])
    U2=np.array([[u[2*nid(i,j)+1] for j in range(nyn)] for i in range(nxn)])
    dx=xs[1]-xs[0]; dy=ys[1]-ys[0]
    def grad(i,j):
        i1=min(i+1,nxn-1); i0=max(i-1,0); j1=min(j+1,nyn-1); j0=max(j-1,0)
        u1x=(U1[i1,j]-U1[i0,j])/(xs[i1]-xs[i0]); u2x=(U2[i1,j]-U2[i0,j])/(xs[i1]-xs[i0])
        u1y=(U1[i,j1]-U1[i,j0])/(ys[j1]-ys[j0]); u2y=(U2[i,j1]-U2[i,j0])/(ys[j1]-ys[j0])
        return u1x,u2x,u1y,u2y,U1[i,j],U2[i,j]
    Ds=E/(1.0-nu**2)
    def WV(i,j):
        u1x,u2x,u1y,u2y,u1v,u2v=grad(i,j)
        e11=u1x;e22=u2y;e12=0.5*(u1y+u2x);tr=e11+e22
        W = h*s0*tr + 0.5*h*Ds*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*tr**2)+0.5*k*(u1v**2+u2v**2)
        s11=s0+Ds*((1-nu)*e11+nu*tr); s22=s0+Ds*((1-nu)*e22+nu*tr)
        s12=Ds*(1-nu)*e12
        return W,(s11,s22,s12),(e11,e22,e12),(u1x,u2x,u1y,u2y)
    iR=np.searchsorted(xs,r); iL=np.searchsorted(xs,-r)
    # J = ∮ (W n1 - sigma_ij n_j u_i,1) ds
    # assume right/left vertical and top/bottom horizontal segments.
    # bottom edge (x2=0, x1 from +r to -r): outward normal of contour
    # (going counterclockwise around tip region inside domain): we select
    # a closed loop; contributions of horizontal edges at x2=0,H for the
    # crack faces add the traction work. To avoid energetic ambiguity we
    # use a contour that does NOT touch the crack faces: rectangle in
    # x1 in [-r, +r], x2 in [eps, H-eps].
    eps=2*dy
    jlo=max(1,int(round(eps/dy))); jhi=nyn-1-jlo
    J=0.0
    # right edge (x1=+r), normal +x1: t1 integral of (W - s11 u1x - s12 u2x) dy? 
    # J_k = ∫ (W n_k - sigma_ij n_j u_i,k) ds
    # k=1: J = ∫ (W n1 - sigma_i1 n1 u_i,1 - sigma_i2 n2 u_i,1) ds
    # on right edge n=(1,0): W - s11 u1x - s12 u2x
    for j in range(jlo,jhi):
        W,s,e,u=WV(iR,j); J += (W - s[0]*u[0] - s[2]*u[1])*dy
    # left edge n=(-1,0): -(W) + s11 u1x + s12 u2x  (with n1=-1)
    for j in range(jlo,jhi):
        W,s,e,u=WV(iL,j); J += (-(W) + s[0]*u[0] + s[2]*u[1])*dy
    # top edge (x2=H-eps) n=(0,1): -(s21 u1x + s22 u2x) ... n2=1
    j=jhi; x1s=range(iL,iR)
    for i in range(iL,iR):
        W,s,e,u=WV(i,j); J += (-(s[2]*u[0] + s[1]*u[1]))*(xs[min(i+1,nxn-1)]-xs[i])
    j=jlo
    for i in range(iL,iR):
        W,s,e,u=WV(i,j); J += ( (s[2]*u[0] + s[1]*u[1]))*(xs[min(i+1,nxn-1)]-xs[i])
    return J

if __name__=="__main__":
    for ell in [1.0]:
        u,xs,ys,nxn,nyn,nid,shape,gauss=solve(ell,400,80,20.0,10.0)
        for r in [0.3,0.5,0.8,1.2]:
            J=Jint(ell,u,xs,ys,nxn,nyn,nid,shape,gauss,r)
            print(f"ell={ell} r={r} J={J:.6f}  J/(s0^2 ell)={J/(s0**2*ell):.6f}")