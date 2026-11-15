import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ---------- material ----------
nu = 1.0/3.0
E1, E2 = 9.0, 1.0
mu1, mu2 = E1/(2*(1+nu)), E2/(2*(1+nu))
kappa = 3-4*nu
Estar1, Estar2 = E1/(1-nu**2), E2/(1-nu**2)
nuprime = nu/(1-nu)

def Dmat(E):
    # plane strain D
    nup = nu
    c = E/((1+nup)*(1-2*nup))
    return c*np.array([[1-nup, nup, 0],[nup,1-nup,0],[0,0,(1-2*nup)/2]])

# Q2 shape functions on [-1,1]^2
def shape(xi,eta):
    Lm=lambda t: 0.5*t*(t-1); L0=lambda t: 1-t*t; Lp=lambda t: 0.5*t*(t+1)
    L = [Lm,Lp,Lp,Lm,L0,Lp,L0,Lm,L0]
    xin=[-1,1,1,-1,0,1,0,-1,0]; etan=[-1,-1,1,1,-1,0,1,0,0]
    N=np.array([L[k](xi if xin[k]!=0 else eta) if False else 0 for k in range(9)],dtype=float)
    # do explicit tensor product
    def vals(t):
        return np.array([0.5*t*(t-1),1-t*t,0.5*t*(t+1)])
    vx=vals(xi); vy=vals(eta)
    # node local coords indices in {-1,0,1}: index map
    ix={-1:0,0:1,1:2}
    N=np.zeros(9)
    for k,(a,b) in enumerate(zip(xin,etan)):
        N[k]=vx[ix[a]]*vy[ix[b]]
    return N

def dshape(xi,eta):
    def vals(t): return np.array([0.5*t*(t-1),1-t*t,0.5*t*(t+1)])
    def dvals(t): return np.array([t-0.5,-2*t,t+0.5])
    vx=vals(xi); vy=vals(eta); dx=dvals(xi); dy=dvals(eta)
    xin=[-1,1,1,-1,0,1,0,-1,0]; etan=[-1,-1,1,1,-1,0,1,0,0]
    ix={-1:0,0:1,1:2}
    Nxi=np.zeros(9); Neta=np.zeros(9)
    for k,(a,b) in enumerate(zip(xin,etan)):
        Nxi[k]=dx[ix[a]]*vy[ix[b]]
        Neta[k]=vx[ix[a]]*dy[ix[b]]
    return Nxi,Neta

# 3x3 Gauss
gp=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)]); gw=np.array([5/9,8/9,5/9])

def geo_grid(L,q,n):
    corners=np.array([0.0]+[ (q**i-1)/(q-1) for i in range(1,n+1)])
    corners*=L/corners[-1]
    mids=np.sqrt(corners[:-1]*corners[1:])
    mids[0]=0.5*corners[1]
    nodes=np.zeros(2*n+1)
    nodes[0::2]=corners; nodes[1::2]=mids
    return nodes

def build_mesh(x0=2e-4, q=1.08, Lx=5.0, H=1.0, h=0.5):
    nx=int(np.ceil(np.log(1+(q-1)*Lx/x0)/np.log(q)))
    xp=geo_grid(Lx,q,nx)  # 0..Lx, length 2nx+1
    x=np.concatenate([-xp[:0:-1], xp])
    def ygrid(L):
        ny=int(np.ceil(np.log(1+(q-1)*L/x0)/np.log(q)))
        y=geo_grid(L,q,ny)
        return y,ny
    y1,ny1=ygrid(h)
    y2,ny2=ygrid(H)
    # nodes for layer 1: x (nxfull) by y1
    # element corners: x corners are indices 0,2,...,2nxfull; y corners 0,2,...
    nxf=len(x); nxcorner=(nxf-1)//2
    ny1c=(len(y1)-1)//2; ny2c=(len(y2)-1)//2
    # global node numbering
    # layer1 nodes: all x grid points, all y1 grid points
    n1=len(x)*len(y1)
    # layer2 nodes: all x grid points, all y2 grid points
    n2=len(x)*len(y2)
    # merge interface nodes: layer1 y index 0, layer2 y index len(y2)-1 (both y=0)
    # for x>=0 merge layer2 top node to layer1 bottom node
    # build a mapping for layer2 node ids -> global
    # We'll assign layer1 ids 0..n1-1; layer2 new ids for unmerged, merged point to layer1 bottom.
    id1=np.arange(n1).reshape(len(y1),len(x))  # [j,i]
    nextid=n1
    id2=np.zeros((len(y2),len(x)),dtype=int)
    # copy non-top rows with new ids
    for j in range(len(y2)):
        for i in range(len(x)):
            if j==len(y2)-1 and x[i]>=0:
                id2[j,i]=id1[0,i]
            else:
                id2[j,i]=nextid; nextid+=1
    ntot=nextid
    # elements
    elems=[]; mat=[]
    # layer1: corner j indices 0,2,..,2ny1c ; i 0,2,..,2nxcorner
    for jj in range(ny1c):
        for ii in range(nxcorner):
            j0=2*jj; i0=2*ii
            conn=[id1[j0,i0],id1[j0,i0+2],id1[j0+2,i0+2],id1[j0+2,i0],
                  id1[j0,i0+1],id1[j0+1,i0+2],id1[j0+2,i0+1],id1[j0+1,i0],
                  id1[j0+1,i0+1]]
            elems.append(conn); mat.append(1)
    # layer2: y2 index 0 corresponds to y=-H, top index len-1 -> y=0
    for jj in range(ny2c):
        for ii in range(nxcorner):
            j0=2*jj; i0=2*ii
            conn=[id2[j0,i0],id2[j0,i0+2],id2[j0+2,i0+2],id2[j0+2,i0],
                  id2[j0,i0+1],id2[j0+1,i0+2],id2[j0+2,i0+1],id2[j0+1,i0],
                  id2[j0+1,i0+1]]
            elems.append(conn); mat.append(2)
    elems=np.array(elems,dtype=int); mat=np.array(mat)
    # coordinates
    coords=np.zeros((ntot,2))
    for j in range(len(y1)):
        for i in range(len(x)):
            coords[id1[j,i]]=[x[i],y1[j]]
    for j in range(len(y2)):
        for i in range(len(x)):
            gid=id2[j,i]
            if j==len(y2)-1 and x[i]>=0:
                continue
            coords[gid]=[x[i],y2[j]-H]  # y2 from 0..H -> -H..0
    return dict(coords=coords,elems=elems,mat=mat,n1=n1,ntot=ntot,
                id1=id1,id2=id2,x=x,y1=y1,y2=y2,nxcorner=nxcorner,ny1c=ny1c,ny2c=ny2c,
                H=H,h=h,Lx=Lx)

def assemble(mesh):
    coords=mesh['coords']; elems=mesh['elems']; mat=mesh['mat']
    ndof=2*mesh['ntot']
    rows=[];cols=[];data=[]
    for e,(conn,m) in enumerate(zip(elems,mat)):
        E=E1 if m==1 else E2
        D=Dmat(E)
        Ke=np.zeros((18,18))
        xy=coords[conn]
        for a in range(3):
            for b in range(3):
                xi=gp[a]; eta=gp[b]; w=gw[a]*gw[b]
                Nxi,Neta=dshape(xi,eta)
                J=np.array([[Nxi@xy[:,0],Nxi@xy[:,1]],[Neta@xy[:,0],Neta@xy[:,1]]])
                detJ=np.linalg.det(J)
                if detJ<=0: raise RuntimeError("bad jac")
                invJ=np.linalg.inv(J)
                dN=np.zeros((2,9))
                dN[0,:]=invJ[0,0]*Nxi+invJ[0,1]*Neta
                dN[1,:]=invJ[1,0]*Nxi+invJ[1,1]*Neta
                B=np.zeros((3,18))
                for k in range(9):
                    B[0,2*k]=dN[0,k]; B[1,2*k+1]=dN[1,k]
                    B[2,2*k]=dN[1,k]; B[2,2*k+1]=dN[0,k]
                Ke+=B.T@D@B*detJ*w
        dofs=np.zeros(18,dtype=int)
        for k in range(9):
            dofs[2*k]=2*conn[k]; dofs[2*k+1]=2*conn[k]+1
        rows.extend(np.repeat(dofs,18)); cols.extend(np.tile(dofs,18)); data.extend(Ke.ravel())
    K=sp.coo_matrix((data,(rows,cols)),shape=(ndof,ndof)).tocsr()
    return K

def apply_loads(mesh,P=1.0,M=0.0):
    coords=mesh['coords']; id1=mesh['id1']; id2=mesh['id2']
    x=mesh['x']; y1=mesh['y1']; y2=mesh['y2']
    H=mesh['H']; h=mesh['h']
    F=np.zeros(2*mesh['ntot'])
    # helper: 1D quadratic edge shape
    def Lvals(t):
        return np.array([0.5*t*(t-1),1-t*t,0.5*t*(t+1)])
    # left boundary i=0 is at x=-Lx
    # layer1: edges along y1, j from 0..len(y1)-2 in steps of 2 (corner intervals)
    for jj in range(mesh['ny1c']):
        j0=2*jj
        # local edge nodes at xi=-1: local 0 (eta=-1),7(eta=0),3(eta=1)
        nodes=[id1[j0,0],id1[j0+1,0],id1[j0+2,0]]
        ys=[y1[j0],y1[j0+1],y1[j0+2]]
        # integrate t_x(y) = -sig11(y)
        for a in range(3):
            xi=gp[a]; w=gw[a]
            L=Lvals(xi)  # shape at eta=xi
            # y and dydeta from the three nodes (quadratic)
            # use Lagrange derivatives
            def Ld(t): return np.array([t-0.5,-2*t,t+0.5])
            dL=Ld(xi)
            yv=sum(L[k]*ys[k] for k in range(3)); dydeta=sum(dL[k]*ys[k] for k in range(3))
            sig=-P/h - 12*M*(yv-h/2)/h**3
            tx=-sig
            for k in range(3):
                F[2*nodes[k]] += L[k]*tx*dydeta*w
    # layer2: y from -H to 0, j0=0..; top at y2=0
    for jj in range(mesh['ny2c']):
        j0=2*jj
        nodes=[id2[j0,0],id2[j0+1,0],id2[j0+2,0]]
        ys=[y2[j0]-H,y2[j0+1]-H,y2[j0+2]-H]  # y coords
        for a in range(3):
            xi=gp[a]; w=gw[a]
            L=Lvals(xi)
            def Ld(t): return np.array([t-0.5,-2*t,t+0.5])
            dL=Ld(xi)
            yv=sum(L[k]*ys[k] for k in range(3)); dydeta=sum(dL[k]*ys[k] for k in range(3))
            # sig11 layer2 = P/H + 12(M+P(h+H)/2)(y+H/2)/H^3 ; P=1,M=0
            sig=P/H + 12*(M + P*(h+H)/2)*(yv+H/2)/H**3
            tx=-sig
            for k in range(3):
                F[2*nodes[k]] += L[k]*tx*dydeta*w
    return F

def solve(mesh,K,F):
    coords=mesh['coords']; Lx=mesh['Lx']; H=mesh['H']; h=mesh['h']
    ndof=2*mesh['ntot']
    # find nodes at right boundary approx Lx, bottom -H and top h
    def find(xx,yy):
        d=np.hypot(coords[:,0]-xx,coords[:,1]-yy)
        return int(np.argmin(d))
    n1=find(Lx,-H); n2=find(Lx,h)
    fixed=[2*n1,2*n1+1,2*n2]
    free=np.setdiff1d(np.arange(ndof),fixed)
    Kff=K[free][:,free]; Kfc=K[free][:,fixed]
    uc=np.zeros(len(fixed))
    rhs=F[free]-Kfc@uc
    u=spla.spsolve(Kff.tocsc(),rhs)
    U=np.zeros(ndof); U[free]=u
    return U
