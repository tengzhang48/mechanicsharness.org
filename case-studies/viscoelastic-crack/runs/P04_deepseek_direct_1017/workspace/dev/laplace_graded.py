import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import time

def make_grid(L,eps,N1=80,N2=240,N3=120):
    t=np.linspace(0,1,N1+1); u=(1-np.cos(np.pi*t))/2
    x1=-L+(L-eps)*u
    t=np.linspace(0,1,N2+1); u=(1-np.cos(np.pi*t))/2
    x2=-eps+eps*u
    t=np.linspace(0,1,N3+1); u=(1-np.cos(np.pi*t))/2
    x3=L*u
    X=np.concatenate([x1,x2[1:],x3[1:]])
    return X

def solve(a,eps,L=12.0,N1=80,N2=240,N3=120,Ny=160,tcoh=1.0):
    X=make_grid(L,eps,N1,N2,N3); Nx=len(X); dY=np.pi/Ny
    def idx(i,j): return (j-1)*(Nx-2)+(i-1)  # j=1..Ny-1
    N=(Ny-1)*(Nx-2)
    rows=[];cols=[];vals=[];b=np.zeros(N)
    # X second derivative coefficients
    ax=np.zeros(Nx); bx=np.zeros(Nx); cx=np.zeros(Nx)
    for i in range(1,Nx-1):
        hm=X[i]-X[i-1]; hp=X[i+1]-X[i]
        ax[i]=2/((hm+hp)*hm); cx[i]=2/((hm+hp)*hp); bx[i]=-(ax[i]+cx[i])
    for j in range(1,Ny):
        for i in range(1,Nx-1):
            r=idx(i,j)
            diag=bx[i]-2/dY**2
            # bottom boundary modification at j=1
            if j==1:
                if X[i]>0:
                    pass  # U0=0, diag -2/dY^2 already
                else:
                    diag+=1/dY**2  # (U2-U1)/dY^2, coefficient -1/dY^2
                    if X[i]>-eps:
                        b[r]+=tcoh/(np.pi*dY)
            rows.append(r);cols.append(r);vals.append(diag)
            if i-1>=1:
                rows.append(r);cols.append(idx(i-1,j));vals.append(ax[i])
            else:
                b[r]-=ax[i]*a
            if i+1<=Nx-2:
                rows.append(r);cols.append(idx(i+1,j));vals.append(cx[i])
            else:
                b[r]-=cx[i]*(a*(j*dY)/np.pi)
            if j-1>=1:
                rows.append(r);cols.append(idx(i,j-1));vals.append(1/dY**2)
            if j+1<=Ny-1:
                rows.append(r);cols.append(idx(i,j+1));vals.append(1/dY**2)
            else:
                b[r]-=(1/dY**2)*a
    A=sp.csr_matrix((vals,(rows,cols)),shape=(N,N))
    U=spla.spsolve(A,b).reshape((Ny-1,Nx-2))
    return X,dY,U

if __name__=='__main__':
    for eps,a in [(0.125,0.4454135),(0.1,0.3993021),(0.15,0.4870175),(0.3,0.679975),(0.5,0.863444)]:
        t0=time.time(); X,dY,U=solve(a,eps,L=12,Ny=160)
        print('eps',eps,'a',a,'time %.2f'%(time.time()-t0),'Nx',len(X))
        # opening at Y=0 for X in cohesive: U0=U1-dY/pi
        # find index near -eps and interpolate
        i=np.argmin(np.abs(X+eps))
        U0=U[0,i-1]-dY/np.pi if X[i]>-eps else U[0,i-1]
        print('  X[idx]',X[i],'U0(-eps approx)',U0,'U at 0-',U[0,np.argmin(np.abs(X-0))-1]-dY/np.pi)
        # sample U at X=-4,-2, y small
        for xq in [-4,-2,-1,-0.5]:
            ii=np.argmin(np.abs(X-xq))
            print('   X',X[ii],'U1',U[0,ii-1],'U0approx',U[0,ii-1]-dY/np.pi)
