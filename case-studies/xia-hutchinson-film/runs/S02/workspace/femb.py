import numpy as np, scipy.sparse as spa, scipy.sparse.linalg as spla

# ======================================================================
# CASE B FEM. Period 2H (set H=1 so period=2). Domain x2 in [-1,1].
# - Existing infinite cracks at x2=+/-1 (period edges): FREE for all x1
#   (outer boundaries of the strip).
# - Additional semi-infinite crack at x2=0 (center), x1<0: free there.
# Far ahead: strip of width 2H, free on outer edges, prestress sigma0 I,
#   plus the central crack only for x1<0.
# Only the additional (central) tip advances.
# Carry prestress as crack-face tractions on ALL free faces (outer edges
# and central crack).
# ======================================================================
nu=0.3;E=1.0;s0=1.0;h=1.0
H=1.0; PER=2*H

def build(ell,nx,ny,Xb,Xf):
    k=E*h/(ell**2*(1.0-nu**2))
    xs=np.linspace(-Xb,Xf,nx+1); ys=np.linspace(-H,H,ny+1)
    nxn=len(xs);nyn=len(ys)
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
                B[0,0::2]=dNxy[0];B[1,1::2]=dNxy[1]
                B[2,0::2]=dNxy[1];B[2,1::2]=dNxy[0]
                ke+=(B.T@D@B)*detJ*h
                Nmat=np.zeros((2,8));Nmat[0,0::2]=N;Nmat[1,1::2]=N
                ke+=(Nmat.T@Nmat)*detJ*k
            for a_ in range(8):
                for b_ in range(8):
                    rows.append(2*nodes[a_//2]+(a_%2));cols.append(2*nodes[b_//2]+(b_%2));vals.append(ke[a_,b_])
    K=spa.csr_matrix((vals,(rows,cols)),shape=(ndof,ndof))
    fixed=[]
    for j in range(nyn):
        n=nid(nx-1,j);fixed+=[2*n,2*n+1]
    # loads: all free faces get t=-sigma0 n.
    # Outer edges x2=+H (j=ny): n=(0,1) -> t=(0,-s0), for ALL x1.
    # Outer edge   x2=-H (j=0): n=(0,-1)-> t=(0,+s0), for ALL x1.
    for i in range(nx):
        Le=xs[i+1]-xs[i]
        F[2*nid(i,ny)+1]+= -s0*Le/2
        F[2*nid(i+1,ny)+1]+= -s0*Le/2
        F[2*nid(i,0)+1]+= s0*Le/2
        F[2*nid(i+1,0)+1]+= s0*Le/2
    return K,F,fixed,xs,ys,nxn,nyn,nid

def solve_tip(ell,nx,ny,Xb,Xf,sh):
    K,F,fixed,xs,ys,nxn,nyn,nid=build(ell,nx,ny,Xb,Xf)
    # central crack x2=0 (j=ny//2): free for x1<0 -> t=(0,?).
    # outward normal from upper half is (0,-1) ? Treat as internal faces:
    # apply equal/opposite tractions on the two coincident faces.
    jm=ny//2
    # redo F for central crack only
    F2=F.copy()
    i0=max(0,int(round((0.0-(-Xb))/(xs[1]-xs[0])))+sh)  # sh>0 -> longer crack
    for i in range(0,i0):
        Le=xs[i+1]-xs[i]
        # lower face (j=jm) belongs to lower half: outward normal (0,+1)?
        # We split x2=0 into two faces. For x1<0 they are free:
        # applied t to cancel prestress. Upper face normal (0,-1)-> t=(0,s0)
        # lower face normal (0,+1) -> t=(0,-s0). Use consistent nodal loads.
        F2[2*nid(i,jm)+1]+= s0*Le/2
        F2[2*nid(i+1,jm)+1]+= s0*Le/2
    ndof=K.shape[0];free=np.ones(ndof,bool);free[fixed]=False
    u=np.zeros(ndof);u[free]=spla.spsolve(K[free][:,free].tocsc(),F2[free])
    # potential energy
    Ds=E/(1.0-nu**2);k=E*h/(ell**2*(1.0-nu**2))
    U1=np.array([[u[2*nid(i,j)] for j in range(nyn)] for i in range(nxn)])
    U2=np.array([[u[2*nid(i,j)+1] for j in range(nyn)] for i in range(nxn)])
    Pi=0.0
    for i in range(nxn):
        for j in range(nyn):
            i1=min(i+1,nxn-1);i0=max(i-1,0);j1=min(j+1,nyn-1);j0=max(j-1,0)
            u1x=(U1[i1,j]-U1[i0,j])/(xs[i1]-xs[i0]);u2x=(U2[i1,j]-U2[i0,j])/(xs[i1]-xs[i0])
            u1y=(U1[i,j1]-U1[i,j0])/(ys[j1]-ys[j0]);u2y=(U2[i,j1]-U2[i,j0])/(ys[j1]-ys[j0])
            e11=u1x;e22=u2y;e12=0.5*(u1y+u2x);tr=e11+e22
            cell=(xs[1]-xs[0])*(ys[1]-ys[0])
            W=h*s0*tr+0.5*h*Ds*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*tr**2)+0.5*k*(U1[i,j]**2+U2[i,j]**2)
            Pi+=W*cell
    return Pi,xs[1]-xs[0]

if __name__=="__main__":
    for r in [0.5,1.0,1.5,2.0]:
        ell=r
        nx,ny,Xb,Xf,sh=600,80,25.0,10.0,10
        P0,dx=solve_tip(ell,nx,ny,Xb,Xf,0)
        P1,dx=solve_tip(ell,nx,ny,Xb,Xf,sh)  # sh>0 tip at -sh*dx => LONGER
        da=sh*dx
        G=-(P1-P0)/da/h
        gam=G*(1/(1-0.3**2))/(s0**2*ell)
        K=np.sqrt(G*E/(s0**2*ell))
        print(f'ell/H={r}: G={G:.6f} gamma_B={gam:.6f} K_B/(s0 sqrt ell)={K:.6f}',flush=True)