import numpy as np, scipy.sparse as spa, scipy.sparse.linalg as spla

nu=0.3;E=1.0;s0=1.0;h=1.0

def build(ell,HH,nx,ny,Xb,Xf):
    k=E*h/(ell**2*(1.0-nu**2))
    xs=np.linspace(-Xb,Xf,nx+1);ys=np.linspace(0.0,HH,ny+1)
    nxn=len(xs);nyn=len(ys)
    nid=lambda i,j:i*nyn+j
    nn=nxn*nyn;ndof=2*nn
    D=(E/(1.0-nu**2))*np.array([[1.0,nu,0.0],[nu,1.0,0.0],[0.0,0.0,(1.0-nu)/2.0]])
    gauss=[(-1/np.sqrt(3),-1/np.sqrt(3)),(1/np.sqrt(3),-1/np.sqrt(3)),
           (1/np.sqrt(3),1/np.sqrt(3)),(-1/np.sqrt(3),1/np.sqrt(3))]
    def shape(xi,eta):
        N=np.array([(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4])
        dN=np.array([[-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4],
                     [-(1-xi)/4,-(1+xi)/4,(1+xi)/4,(1-xi)/4]])
        return N,dN
    rows=[];cols=[];vals=[]
    for i in range(nx):
        for j in range(ny):
            nodes=[nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)]
            xy=np.array([[xs[i],ys[j]],[xs[i+1],ys[j]],[xs[i+1],ys[j+1]],[xs[i],ys[j+1]]])
            ke=np.zeros((8,8))
            for xi,eta in gauss:
                N,dN=shape(xi,eta)
                J=dN@xy;detJ=np.linalg.det(J)
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
    return K,xs,ys,nxn,nyn,nid

def solve(ell,HH,nx,ny,Xb,Xf,sh):
    K,xs,ys,nxn,nyn,nid=build(ell,HH,nx,ny,Xb,Xf)
    ndof=K.shape[0];F=np.zeros(ndof)
    fixed=[]
    for j in range(nyn):
        n=nid(nx-1,j);fixed+=[2*n,2*n+1]
    dx=xs[1]-xs[0]
    itip=int(round((0.0-(-Xb))/dx))+sh
    for i in range(0,itip):
        Le=xs[i+1]-xs[i]
        F[2*nid(i,0)+1]+= s0*Le/2
        F[2*nid(i+1,0)+1]+= s0*Le/2
        F[2*nid(i,ny)+1]+= -s0*Le/2
        F[2*nid(i+1,ny)+1]+= -s0*Le/2
    free=np.ones(ndof,bool);free[fixed]=False
    u=np.zeros(ndof);u[free]=spla.spsolve(K[free][:,free].tocsc(),F[free])
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
            cell=dx*(ys[1]-ys[0])
            W=h*s0*tr+0.5*h*Ds*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*tr**2)+0.5*k*(U1[i,j]**2+U2[i,j]**2)
            Pi+=W*cell
    return Pi,dx

if __name__=="__main__":
    # Isolated crack limit: ell=1 fixed, H large.
    ell=1.0
    for HH in [1.0,2.0,4.0,8.0]:
        nx,ny,sh=700,80,10
        P0,dx=solve(ell,HH,nx,ny,25.0,10.0,0)
        P1,_=solve(ell,HH,nx,ny,25.0,10.0,sh)
        G=-(P1-P0)/(sh*dx)/h
        print(f'H/ell={HH/ell:.1f}: G_A={G:.6f} gamma_A=G*Ebar/(s0^2 ell)={G*(1/(1-0.3**2)):.6f} K_A={np.sqrt(G):.6f}',flush=True)