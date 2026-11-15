import numpy as np, fem

def jdomain(coords,elems,mat,U,r1=0.05,r2=0.2):
    total_pos=0.0; total_neg=0.0
    gp=fem.gp; gw=fem.gw
    for e,(conn,m) in enumerate(zip(elems,mat)):
        xy=coords[conn]
        E=fem.E1 if m==1 else fem.E2
        D=fem.Dmat(E)
        ue=np.zeros(18)
        for k in range(9):
            ue[2*k]=U[2*conn[k]]; ue[2*k+1]=U[2*conn[k]+1]
        for a in range(3):
            for b in range(3):
                xi=gp[a]; eta=gp[b]; w=gw[a]*gw[b]
                N=fem.shape(xi,eta); Nxi,Neta=fem.dshape(xi,eta)
                J=np.array([[Nxi@xy[:,0],Nxi@xy[:,1]],[Neta@xy[:,0],Neta@xy[:,1]]])
                detJ=np.linalg.det(J); invJ=np.linalg.inv(J)
                dN=np.zeros((2,9))
                dN[0,:]=invJ[0,0]*Nxi+invJ[0,1]*Neta
                dN[1,:]=invJ[1,0]*Nxi+invJ[1,1]*Neta
                eps=np.zeros(3); ux1=0.0; uy1=0.0
                for k in range(9):
                    eps[0]+=dN[0,k]*ue[2*k]
                    eps[1]+=dN[1,k]*ue[2*k+1]
                    eps[2]+=dN[1,k]*ue[2*k]+dN[0,k]*ue[2*k+1]
                    ux1+=dN[0,k]*ue[2*k]; uy1+=dN[0,k]*ue[2*k+1]
                sig=D@eps
                W=0.5*(sig[0]*eps[0]+sig[1]*eps[1]+2*sig[2]*eps[2])
                x=N@xy[:,0]; y=N@xy[:,1]; r=np.hypot(x,y)
                if r<=r1: q=1.0; dq=0.0
                elif r>=r2: continue
                else:
                    t=(r-r1)/(r2-r1); q=1-3*t*t+2*t**3; dq=(-6*t+6*t*t)/(r2-r1)
                qx=dq*x/r; qy=dq*y/r
                # sigma_ij u_i,1 : (sig_xx ux1 + sig_xy uy1, sig_xy ux1+sig_yy uy1) dot qgrad
                term=(sig[0]*ux1+sig[2]*uy1)*qx+(sig[2]*ux1+sig[1]*uy1)*qy
                wj=term - W*qx
                total_pos+=wj*detJ*w
                total_neg+=(W*qx-term)*detJ*w
    return total_pos,total_neg

if __name__=='__main__':
    d=np.load('fine_data.npz')
    coords,elems,mat,U=d['coords'],d['elems'],d['mat'],d['U']
    for r1,r2 in [(0.02,0.1),(0.05,0.2),(0.1,0.3),(0.05,0.5)]:
        p,n=jdomain(coords,elems,mat,U,r1,r2)
        print(r1,r2,"pos",p,"neg",n)
