import numpy as np
import fem
exec(open('analytic.py').read().split('# ---------------- analytic Cartesian')[0])

def airy_stress_pt(mat,x,y):
    c=c1 if mat==1 else c2
    r=np.hypot(x,y); th=np.arctan2(y,x)
    f=c@basis(th,s); fp=c@dbasis(th,s); fpp=c@ddbasis(th,s)
    Prr=(s+2)*f+fpp; Ptt=(s+2)*(s+1)*f; Prt=-(s+1)*fp
    rs=r**s
    Sxx=rs*(Prr*np.cos(th)**2+Ptt*np.sin(th)**2-2*Prt*np.sin(th)*np.cos(th))
    Syy=rs*(Prr*np.sin(th)**2+Ptt*np.cos(th)**2+2*Prt*np.sin(th)*np.cos(th))
    Sxy=rs*((Prr-Ptt)*np.sin(th)*np.cos(th)+Prt*(np.cos(th)**2-np.sin(th)**2))
    return np.real(Sxx),np.real(Syy),np.real(Sxy)

def edge_load(mesh,F,conn_edge,normal):
    """conn_edge: 3 global node ids at eta=-1,0,+1 along edge;
       normal: (nx,ny) outward. Integrate analytic traction."""
    coords=mesh['coords']
    xy=coords[conn_edge]
    # material from midpoint y: if >0 layer1 else layer2 (edges are within one layer)
    ymid=xy[1,1]
    mat=1 if ymid>1e-9 else 2
    def Lvals(t): return np.array([0.5*t*(t-1),1-t*t,0.5*t*(t+1)])
    def Lder(t): return np.array([t-0.5,-2*t,t+0.5])
    for a in range(3):
        xi=fem.gp[a]; w=fem.gw[a]
        L=Lvals(xi); dL=Lder(xi)
        x=sum(L[k]*xy[k,0] for k in range(3)); y=sum(L[k]*xy[k,1] for k in range(3))
        dxd=np.array([dL[0]*xy[0,0]+dL[1]*xy[1,0]+dL[2]*xy[2,0],
                      dL[0]*xy[0,1]+dL[1]*xy[1,1]+dL[2]*xy[2,1]])
        ds=np.linalg.norm(dxd)
        # ensure positive orientation (not important, integrate abs)
        sxx,syy,sxy=airy_stress_pt(mat,x,y)
        tx=normal[0]*sxx+normal[1]*sxy
        ty=normal[0]*sxy+normal[1]*syy
        for k in range(3):
            F[2*conn_edge[k]] += L[k]*tx*ds*w
            F[2*conn_edge[k]+1] += L[k]*ty*ds*w

def apply_tractions(mesh):
    F=np.zeros(2*mesh['ntot'])
    id1=mesh['id1']; id2=mesh['id2']; nx=mesh['nxcorner']; ny1=mesh['ny1c']; ny2=mesh['ny2c']
    # top layer1: j0=2*ny1-2, edge eta=+1 nodes [3,6,2]
    for ii in range(nx):
        i0=2*ii; j0=2*ny1-2
        conn=[id1[j0+2,i0],id1[j0+2,i0+1],id1[j0+2,i0+2]]
        edge_load(mesh,F,conn,(0,1))
    # bottom layer2: elements j0=0, edge eta=-1 nodes [0,4,1]
    for ii in range(nx):
        i0=2*ii
        conn=[id2[0,i0],id2[0,i0+1],id2[0,i0+2]]
        edge_load(mesh,F,conn,(0,-1))
    # left/right edges: layer1 and layer2, elements along y
    for layer in [1,2]:
        ny=ny1 if layer==1 else ny2
        idg=id1 if layer==1 else id2
        for jj in range(ny):
            j0=2*jj
            # left edge xi=-1 nodes [0,7,3] at i0=0
            conn=[idg[j0,0],idg[j0+1,0],idg[j0+2,0]]
            edge_load(mesh,F,conn,(-1,0))
            # right edge xi=+1 nodes [1,5,2] at i0=2*nx-2
            i0=2*nx-2
            conn=[idg[j0,i0+2],idg[j0+1,i0+2],idg[j0+2,i0+2]]
            edge_load(mesh,F,conn,(1,0))
    return F

def run(x0=0.005,q=1.2,a=1.0):
    m=fem.build_mesh(x0=x0,q=q,Lx=a,H=a,h=a)
    K=fem.assemble(m); F=apply_tractions(m)
    print("net force",F[0::2].sum(),F[1::2].sum())
    U=fem.solve(m,K,F)
    return m,U
