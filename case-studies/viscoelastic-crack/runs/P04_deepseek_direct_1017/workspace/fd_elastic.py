import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def solve_elastic(lc, Delta, G=1.0, tauc=1.0, L=8.0, Nx=1200, Ny=240, h=1.0):
    # grid in x: cluster near 0 and -lc using a tanh mapping combined
    # simple: uniform grid but with extra points near 0 and -lc via a piecewise map
    # use a global coordinate s in [-1,1], x = L * f(s) with f odd; cluster at 0
    s = np.linspace(-1,1,Nx)
    # map s->x with clustering at 0: x = L * sign(s)*|s|^2? gives clustering at 0
    x = L*np.sign(s)*np.abs(s)**2
    # also cluster near -lc by adding a node? skip for now
    y = h*(1-np.cos(np.pi*np.arange(Ny)/ (Ny-1)))/2  # cluster at both ends
    X,Y = np.meshgrid(x,y,indexing='ij')
    Nx,Ny = X.shape
    idx = lambda i,j: i*Ny+j
    rows=[]; cols=[]; vals=[]; b=np.zeros(Nx*Ny)
    dx = np.diff(x); dy=np.diff(y)
    # use finite volume-ish coefficients
    for i in range(1,Nx-1):
        hxL=x[i]-x[i-1]; hxR=x[i+1]-x[i]
        for j in range(1,Ny-1):
            hyL=y[j]-y[j-1]; hyR=y[j+1]-y[j]
            cL= 2/(hxL*(hxL+hxR)); cR= 2/(hxR*(hxL+hxR))
            cD= 2/(hyL*(hyL+hyR)); cU= 2/(hyR*(hyL+hyR))
            c0= -(cL+cR+cD+cU)
            k=idx(i,j)
            rows += [k,k,k,k,k]; cols += [idx(i-1,j),idx(i+1,j),idx(i,j-1),idx(i,j+1),k]
            vals += [cL,cR,cD,cU,c0]
            b[k]=0.0
    # boundaries x = +-L
    for j in range(Ny):
        k=idx(0,j); rows.append(k); cols.append(k); vals.append(1.0)
        b[k]=Delta if x[0]<0 else Delta*y[j]/h  # left far behind: w=Delta
        k=idx(Nx-1,j); rows.append(k); cols.append(k); vals.append(1.0)
        b[k]=Delta*y[j]/h
    # y=h top
    for i in range(Nx):
        k=idx(i,Ny-1); rows.append(k); cols.append(k); vals.append(1.0); b[k]=Delta
    # y=0 interface
    for i in range(Nx):
        k=idx(i,0)
        if x[i]>0:
            rows.append(k); cols.append(k); vals.append(1.0); b[k]=0.0
        else:
            # Neumann: (w[i,1]-w[i,0])/dy0 = q
            dy0=y[1]-y[0]
            if x[i] > -lc:
                q=tauc/G
            else:
                q=0.0
            rows += [k,k]; cols += [idx(i,1),k]; vals += [1.0/dy0,-1.0/dy0]; b[k]=q
    A=sp.csr_matrix((vals,(rows,cols)),shape=(Nx*Ny,Nx*Ny))
    w=spla.spsolve(A,b).reshape(Nx,Ny)
    # stress tau2 at y=0+ for x>0: G*(w[:,1]-w[:,0])/dy0 (but w[:,0]=0 there)
    tau2=G*(w[:,1]-w[:,0])/(y[1]-y[0])
    # opening at -lc
    U=np.interp(-lc,x,w[:,0])
    return x,y,w,tau2,U

if __name__=='__main__':
    for lc in [0.1,0.2,0.3,0.5]:
        x,y,w,tau2,U=solve_elastic(lc,0.5,G=1,tauc=1)
        # value near x=0+
        i0=np.searchsorted(x,0)
        print('lc',lc,'tau2(0+)',tau2[i0:i0+3],'U(-lc)',U,'Delta',0.5,'a',0.5)
