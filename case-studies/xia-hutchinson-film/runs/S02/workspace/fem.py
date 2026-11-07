import numpy as np
import scipy.sparse as spa
import scipy.sparse.linalg as spla

# ======================================================================
# CASE A -- bilinear-quad FEM of one period strip [0,H], plane stress +
# elastic foundation, prestress applied as crack-face tractions.
# Extract near-tip opening and K, then G = K^2/E.
# ======================================================================
nu   = 0.3
E    = 1.0
s0   = 1.0
h    = 1.0
H    = 1.0

def solve_caseA(ell, nx=600, ny=60, Xback=25.0, Xfront=8.0):
    k = E*h/(ell**2*(1.0-nu**2))
    xs = np.linspace(-Xback, Xfront, nx+1)
    ys = np.linspace(0.0, H, ny+1)
    nxn=len(xs); nyn=len(ys)
    def nid(i,j): return i*nyn+j
    nn=nxn*nyn
    ndof=2*nn
    D = (E/(1.0-nu**2))*np.array([[1.0,nu,0.0],[nu,1.0,0.0],[0.0,0.0,(1.0-nu)/2.0]])
    gauss=[(-1/np.sqrt(3),-1/np.sqrt(3)),(1/np.sqrt(3),-1/np.sqrt(3)),
           (1/np.sqrt(3),1/np.sqrt(3)),(-1/np.sqrt(3),1/np.sqrt(3))]
    def shape(xi,eta):
        N=np.array([(1-xi)*(1-eta)/4,(1+xi)*(1-eta)/4,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4])
        dN=np.array([[-(1-eta)/4,(1-eta)/4,(1+eta)/4,-(1+eta)/4],
                     [-(1-xi)/4,-(1+xi)/4,(1+xi)/4,(1-xi)/4]])
        return N,dN
    rows=[];cols=[];vals=[];F=np.zeros(ndof)
    def add(a,b,v): rows.append(a);cols.append(b);vals.append(v)
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
            for a in range(8):
                for b in range(8):
                    add(2*nodes[a//2]+(a%2),2*nodes[b//2]+(b%2),ke[a,b])
    K=spa.csr_matrix((vals,(rows,cols)),shape=(ndof,ndof))
    # Dirichlet far ahead
    fixed=[]
    for j in range(nyn):
        n=nid(nx-1,j); fixed+=[2*n,2*n+1]
    # crack-face loads: x2=0, x1<0 -> t=(0,s0); x2=H, x1<0 -> t=(0,-s0)
    i0=np.searchsorted(xs,0.0)
    for i in range(0,i0):
        Le=xs[i+1]-xs[i]
        # bottom face
        F[2*nid(i,0)+1]+= s0*Le/2
        F[2*nid(i+1,0)+1]+= s0*Le/2
        # top face (outward normal +y, prestress traction +s0 y -> apply -s0)
        F[2*nid(i,ny)+1]+= -s0*Le/2
        F[2*nid(i+1,ny)+1]+= -s0*Le/2
    return K,F,fixed,nxn,nyn,xs,ys,nid

def run(ell):
    K,F,fixed,nxn,nyn,xs,ys,nid=solve_caseA(ell)
    ndof=K.shape[0]
    free=np.ones(ndof,bool); free[fixed]=False
    Kf=K[free][:,free]; Ff=F[free]
    u=np.zeros(ndof)
    u[free]=spla.spsolve(Kf.tocsc(),Ff)
    return u,nxn,nyn,xs,ys,nid

if __name__=="__main__":
    ell=1.0
    u,nxn,nyn,xs,ys,nid=run(ell)
    # opening at crack face x2=0: u2 (bottom) - u2 (top), for x1<0
    print("xs range",xs[0],xs[-1])
    # compute opening just behind tip for decreasing x1<0
    i0=np.searchsorted(xs,0.0)
    openb=[]
    for i in range(i0-1,i0-8,-1):
        b=nid(i,0); t=nid(i,nyn-1)
        openb.append((xs[i], u[2*t+1]-u[2*b+1]))
    for x,o in openb: print(f"x={x: .4f} open={o: .6e}")
    print("max u2 bottom near tip", max(abs(u[2*nid(i,0)+1]) for i in range(i0-8,i0)))