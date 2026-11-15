import numpy as np
import fem

def elem_stress(mesh,U,e,xi,eta):
    conn=mesh['elems'][e]; m=mesh['mat'][e]
    xy=mesh['coords'][conn]
    Nxi,Neta=fem.dshape(xi,eta)
    J=np.array([[Nxi@xy[:,0],Nxi@xy[:,1]],[Neta@xy[:,0],Neta@xy[:,1]]])
    detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
    dN=np.zeros((2,9))
    dN[0,:]=invJ[0,0]*Nxi+invJ[0,1]*Neta
    dN[1,:]=invJ[1,0]*Nxi+invJ[1,1]*Neta
    ue=np.zeros(18)
    for k in range(9):
        ue[2*k]=U[2*conn[k]]; ue[2*k+1]=U[2*conn[k]+1]
    eps=np.zeros(3)
    for k in range(9):
        eps[0]+=dN[0,k]*ue[2*k]
        eps[1]+=dN[1,k]*ue[2*k+1]
        eps[2]+=dN[1,k]*ue[2*k]+dN[0,k]*ue[2*k+1]
    E=fem.E1 if m==1 else fem.E2
    sig=fem.Dmat(E)@eps
    # coordinates
    N=fem.shape(xi,eta); x=N@xy[:,0]; y=N@xy[:,1]
    return x,y,sig[0],sig[1],sig[2]

def interface_data(mesh,U):
    """return arrays r, sig22_1, sig12_1, sig22_2, sig12_2 for interface x>0"""
    elems=mesh['elems']; mat=mesh['mat']
    xi_g=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)])
    out=[]
    # layer1 elements adjacent to interface: need identify by node coords
    for e,(conn,m) in enumerate(zip(elems,mat)):
        xy=mesh['coords'][conn]
        # check bottom edge midpoint y ~ 0 and material 1, and x>0 somewhere
        ybot=0.5*(xy[0,1]+xy[1,1])
        if m==1 and abs(ybot)<1e-9:
            for xi in xi_g:
                x,y,sxx,syy,sxy=elem_stress(mesh,U,e,xi,-1.0)
                out.append((x,m,syy,sxy))
        # layer2 top edge: local nodes 3,2 have y=0 (eta=+1)
        ytop=0.5*(xy[2,1]+xy[3,1])
        if m==2 and abs(ytop)<1e-9:
            for xi in xi_g:
                x,y,sxx,syy,sxy=elem_stress(mesh,U,e,xi,1.0)
                out.append((x,m,syy,sxy))
    out=np.array(out)
    r1=out[(out[:,1]==1)&(out[:,0]>0)]
    r2=out[(out[:,1]==2)&(out[:,0]>0)]
    # sort and match by r
    r1=r1[np.argsort(r1[:,0])]; r2=r2[np.argsort(r2[:,0])]
    return r1,r2

def crack_jump(mesh,U):
    """displacement jump at crack faces x<0; pair layer1 bottom edge and layer2 top edge at same x"""
    elems=mesh['elems']; mat=mesh['mat']; coords=mesh['coords']
    xi_g=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)])
    d1={}; d2={}
    for e,(conn,m) in enumerate(zip(elems,mat)):
        xy=coords[conn]
        if m==1:
            ybot=0.5*(xy[0,1]+xy[1,1])
            if abs(ybot)<1e-9:
                for xi in xi_g:
                    N=fem.shape(xi,-1.0); x=N@xy[:,0]
                    if x<0:
                        ux=sum(N[k]*U[2*conn[k]] for k in range(9))
                        uy=sum(N[k]*U[2*conn[k]+1] for k in range(9))
                        d1[round(x,12)]=(ux,uy)
        else:
            ytop=0.5*(xy[2,1]+xy[3,1])
            if abs(ytop)<1e-9:
                for xi in xi_g:
                    N=fem.shape(xi,1.0); x=N@xy[:,0]
                    if x<0:
                        ux=sum(N[k]*U[2*conn[k]] for k in range(9))
                        uy=sum(N[k]*U[2*conn[k]+1] for k in range(9))
                        d2[round(x,12)]=(ux,uy)
    keys=sorted(set(d1)&set(d2))
    r=np.array([-k for k in keys])
    dx=np.array([d1[k][0]-d2[k][0] for k in keys])
    dy=np.array([d1[k][1]-d2[k][1] for k in keys])
    return r,dx,dy
