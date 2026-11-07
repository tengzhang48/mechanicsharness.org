"""
Independent 2D J-integral verification, generalized to Case A and Case B.

Half-strip symmetry reduction.  Domain x2 in [0,Wd], x1 in [-Lx,Lx].
Horizontal-edge boundary types:
  'mix'  : x1<i0 traction-free crack, x1>=i0 symmetry (u2=0)
  'free' : traction-free crack for all x1
Far fields are Dirichlet from the stated bounded equilibria.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import quad

nu=0.3; Ebar=1.0; s0=1.0; ell=1.0; h=1.0; k=1.0

def solve(Wd, Lx, Nx, Ny, left_u2, right_u2, bc0, bcN):
    """Return U (shape (Nx+1,Ny+1,2)), dx, dy."""
    dx=2*Lx/Nx; dy=Wd/Ny; i0=Nx//2
    c1=(1-nu)/2.0; c2=(1+nu)/2.0; ell2=ell**2
    nn=(Nx+1)*(Ny+1); dof=2*nn
    def nid(i,j): return i*(Ny+1)+j
    rows=[]; cols=[]; vals=[]; b=np.zeros(dof)
    def add(r,c,v): rows.append(r); cols.append(c); vals.append(v)
    for i in range(Nx+1):
        for j in range(Ny+1):
            n=nid(i,j); r1=2*n; r2=2*n+1
            bot=(j==0); top=(j==Ny)
            if i==0:
                add(r1,r1,1.0); add(r2,r2,1.0); b[r1]=0.0; b[r2]=left_u2[j]
            elif i==Nx:
                add(r1,r1,1.0); add(r2,r2,1.0); b[r1]=0.0; b[r2]=right_u2[j]
            elif bot or top:
                btype = bc0 if bot else bcN
                is_free = (btype=='free') or (btype=='mix' and i<i0)
                sgn = +1.0 if bot else -1.0
                jn = 1 if bot else Ny-1
                if is_free:
                    add(r1,2*nid(i,jn),sgn/dy); add(r1,2*n,-sgn/dy)
                    add(r1,2*nid(i+1,j)+1,1.0/(2*dx)); add(r1,2*nid(i-1,j)+1,-1.0/(2*dx))
                    add(r2,2*nid(i,jn)+1,sgn/dy); add(r2,2*n+1,-sgn/dy)
                    add(r2,2*nid(i+1,j),nu/(2*dx)); add(r2,2*nid(i-1,j),-nu/(2*dx))
                    b[r2]=-s0/Ebar
                else:
                    add(r2,r2,1.0); b[r2]=0.0
                    add(r1,2*nid(i,jn),sgn/dy); add(r1,2*n,-sgn/dy)
                    add(r1,2*nid(i+1,j)+1,1.0/(2*dx)); add(r1,2*nid(i-1,j)+1,-1.0/(2*dx))
            else:
                ip,im=nid(i+1,j),nid(i-1,j); jp,jm=nid(i,j+1),nid(i,j-1)
                add(r1,r1,-2/dx**2-(1-nu)/dy**2-1/ell2)
                add(r1,2*ip,1/dx**2); add(r1,2*im,1/dx**2)
                add(r1,2*jp,c1/dy**2); add(r1,2*jm,c1/dy**2)
                add(r1,2*nid(i+1,j+1)+1, c2/(4*dx*dy)); add(r1,2*nid(i+1,j-1)+1,-c2/(4*dx*dy))
                add(r1,2*nid(i-1,j+1)+1,-c2/(4*dx*dy)); add(r1,2*nid(i-1,j-1)+1, c2/(4*dx*dy))
                add(r2,r2,-c1*2/dx**2-2/dy**2-1/ell2)
                add(r2,2*ip+1,c1/dx**2); add(r2,2*im+1,c1/dx**2)
                add(r2,2*jp+1,1/dy**2); add(r2,2*jm+1,1/dy**2)
                add(r2,2*nid(i+1,j+1), c2/(4*dx*dy)); add(r2,2*nid(i+1,j-1),-c2/(4*dx*dy))
                add(r2,2*nid(i-1,j+1),-c2/(4*dx*dy)); add(r2,2*nid(i-1,j-1), c2/(4*dx*dy))
    A=sp.csr_matrix((vals,(rows,cols)),shape=(dof,dof))
    sol=spla.spsolve(A,b)
    return sol.reshape(Nx+1,Ny+1,2), dx, dy

def Jarc(U, dx, dy, r0, x2c, ysign=+1):
    """J over a semicircle of radius r0 around (x1=0, x2=x2c), on the side
    x2>x2c (ysign=+1).  Domain fields in U (Nx+1,Ny+1,2), x1=-Lx..Lx,
    x2=0..Wd.  Returns the arc J (energy/length)."""
    Nx=U.shape[0]-1; Ny=U.shape[1]-1; Lx=Nx*dx/2.0
    u1=U[:,:,0]; u2=U[:,:,1]
    def ddx(f):
        g=np.zeros_like(f); g[1:-1,:]=(f[2:,:]-f[:-2,:])/(2*dx)
        g[0,:]=(f[1,:]-f[0,:])/dx; g[-1,:]=(f[-1,:]-f[-2,:])/dx; return g
    def ddy(f):
        g=np.zeros_like(f); g[:,1:-1]=(f[:,2:]-f[:,:-2])/(2*dy)
        g[:,0]=(f[:,1]-f[:,0])/dy; g[:,-1]=(f[:,-1]-f[:,-2])/dy; return g
    u1x,u1y=ddx(u1),ddy(u1); u2x,u2y=ddx(u2),ddy(u2)
    def interp(field,xp,yp):
        fi=(xp+Lx)/dx; fj=(yp)/dy
        i=int(np.floor(fi)); j=int(np.floor(fj))
        i=min(max(i,0),Nx-1); j=min(max(j,0),Ny-1)
        tx=fi-i; ty=fj-j
        return ((1-tx)*(1-ty)*field[i,j]+tx*(1-ty)*field[i+1,j]
                +(1-tx)*ty*field[i,j+1]+tx*ty*field[i+1,j+1])
    def integrand(th):
        xp=r0*np.cos(th); yp=x2c+ysign*r0*np.sin(th)
        e11=interp(u1x,xp,yp); e22=interp(u2y,xp,yp)
        e12=0.5*(interp(u1y,xp,yp)+interp(u2x,xp,yp)); egg=e11+e22
        w=h*s0*egg+0.5*h*Ebar*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*egg**2)
        w+=0.5*k*(interp(u1,xp,yp)**2+interp(u2,xp,yp)**2)
        s11=s0+Ebar*((1-nu)*e11+nu*egg); s22=s0+Ebar*((1-nu)*e22+nu*egg)
        s12=Ebar*(1-nu)*e12
        nx=np.cos(th); ny=ysign*np.sin(th)
        T1=h*(s11*nx+s12*ny); T2=h*(s12*nx+s22*ny)
        d1=interp(u1x,xp,yp); d2=interp(u2x,xp,yp)
        return (w*nx-(T1*d1+T2*d2))*r0
    val,_=quad(integrand,0.0,np.pi,limit=400)
    return val
