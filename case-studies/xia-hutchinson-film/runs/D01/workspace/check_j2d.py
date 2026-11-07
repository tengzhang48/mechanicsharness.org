"""
Independent 2D verification of the analytic energy-release results.

Solves the plane-stress film-on-foundation system in the symmetry half-strip
x2 in [0,H], with a semi-infinite crack on x2=0 (x1<0), using finite
differences with duplicated crack-face treatment (traction-free) and a
symmetry plane on x2=0 (x1>0).  Then evaluates the J-integral on a circular
arc around the tip; J_full = 2*J_arc (symmetry) and G = J_full/h.

This route does NOT use the analytic wake energy-difference formula.

Bounded by `timeout`.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import quad

# ---- parameters (nondimensional: sigma0=Ebar=k=h=ell=1, nu=0.3) ----
nu = 0.3
Ebar = 1.0
s0 = 1.0
ell = 1.0
h = 1.0
k = 1.0

H = 1.0
Lx = 4.0
Nx = 400
Ny = 50
dx = 2*Lx/Nx
dy = H/Ny
i0 = Nx//2                     # x1=0
x1 = -Lx + np.arange(Nx+1)*dx
x2 = np.arange(Ny+1)*dy

c1 = (1-nu)/2.0
c2 = (1+nu)/2.0
ell2 = ell**2

nn = (Nx+1)*(Ny+1)
dof = 2*nn
def nid(i,j):
    return i*(Ny+1)+j

rows=[]; cols=[]; vals=[]; b=np.zeros(dof)
def add(r,c,v):
    rows.append(r); cols.append(c); vals.append(v)

Bwake = -s0*ell/(Ebar*np.cosh(H/(2*ell)))
wake_u2 = Bwake*np.sinh((x2 - H/2)/ell)

for i in range(Nx+1):
    for j in range(Ny+1):
        n = nid(i,j); r1 = 2*n; r2 = 2*n+1
        if i==0:                                   # far behind: wake Dirichlet
            add(r1,r1,1.0); add(r2,r2,1.0)
            b[r1]=0.0; b[r2]=wake_u2[j]
        elif i==Nx:                                # far ahead: u=0
            add(r1,r1,1.0); add(r2,r2,1.0)
            b[r1]=0.0; b[r2]=0.0
        elif j==0:
            if i<i0:                               # crack face (traction free)
                add(r1, 2*nid(i,1), 1.0/dy); add(r1, 2*nid(i,0), -1.0/dy)
                add(r1, 2*nid(i+1,0)+1, 1.0/(2*dx)); add(r1, 2*nid(i-1,0)+1, -1.0/(2*dx))
                add(r2, 2*nid(i,1)+1, 1.0/dy); add(r2, 2*nid(i,0)+1, -1.0/dy)
                add(r2, 2*nid(i+1,0), nu/(2*dx)); add(r2, 2*nid(i-1,0), -nu/(2*dx))
                b[r2] = -s0/Ebar
            else:                                  # symmetry plane u2=0
                add(r2,r2,1.0); b[r2]=0.0
                add(r1, 2*nid(i,1), 1.0/dy); add(r1, 2*nid(i,0), -1.0/dy)
                add(r1, 2*nid(i+1,0)+1, 1.0/(2*dx)); add(r1, 2*nid(i-1,0)+1, -1.0/(2*dx))
        elif j==Ny:
            if i<i0:
                add(r1, 2*nid(i,Ny), 1.0/dy); add(r1, 2*nid(i,Ny-1), -1.0/dy)
                add(r1, 2*nid(i+1,Ny)+1, 1.0/(2*dx)); add(r1, 2*nid(i-1,Ny)+1, -1.0/(2*dx))
                add(r2, 2*nid(i,Ny)+1, 1.0/dy); add(r2, 2*nid(i,Ny-1)+1, -1.0/dy)
                add(r2, 2*nid(i+1,Ny), nu/(2*dx)); add(r2, 2*nid(i-1,Ny), -nu/(2*dx))
                b[r2] = -s0/Ebar
            else:
                add(r2,r2,1.0); b[r2]=0.0
                add(r1, 2*nid(i,Ny), 1.0/dy); add(r1, 2*nid(i,Ny-1), -1.0/dy)
                add(r1, 2*nid(i+1,Ny)+1, 1.0/(2*dx)); add(r1, 2*nid(i-1,Ny)+1, -1.0/(2*dx))
        else:
            ip, im = nid(i+1,j), nid(i-1,j)
            jp, jm = nid(i,j+1), nid(i,j-1)
            ipp,jpp = 2*nid(i+1,j+1), 2*nid(i+1,j+1)
            # component 1
            add(r1, r1, -2/dx**2 - (1-nu)/dy**2 - 1/ell2)
            add(r1, 2*ip, 1/dx**2); add(r1, 2*im, 1/dx**2)
            add(r1, 2*jp, c1/dy**2); add(r1, 2*jm, c1/dy**2)
            add(r1, 2*nid(i+1,j+1)+1,  c2/(4*dx*dy))
            add(r1, 2*nid(i+1,j-1)+1, -c2/(4*dx*dy))
            add(r1, 2*nid(i-1,j+1)+1, -c2/(4*dx*dy))
            add(r1, 2*nid(i-1,j-1)+1,  c2/(4*dx*dy))
            # component 2
            add(r2, r2, -c1*2/dx**2 - 2/dy**2 - 1/ell2)
            add(r2, 2*ip+1, c1/dx**2); add(r2, 2*im+1, c1/dx**2)
            add(r2, 2*jp+1, 1/dy**2); add(r2, 2*jm+1, 1/dy**2)
            add(r2, 2*nid(i+1,j+1),  c2/(4*dx*dy))
            add(r2, 2*nid(i+1,j-1), -c2/(4*dx*dy))
            add(r2, 2*nid(i-1,j+1), -c2/(4*dx*dy))
            add(r2, 2*nid(i-1,j-1),  c2/(4*dx*dy))

A = sp.csr_matrix((vals,(rows,cols)), shape=(dof,dof))
print("assembled", A.shape, "nnz", A.nnz)
sol = spla.spsolve(A, b)
U = sol.reshape(nn,2)              # U[node, comp]

u1 = U[:,0].reshape(Nx+1,Ny+1)
u2 = U[:,1].reshape(Nx+1,Ny+1)

# ---- solver sanity: compare to analytic wake at x1 near left boundary ----
iL = 2
err_wake = np.max(np.abs(u2[iL,:]-wake_u2))
print(f"solver check: max|u2_num - u2_wake| at i={iL}: {err_wake:.3e}")
print(f"solver check: max|u| in far-ahead column i=Nx-2: {np.max(np.abs(U.reshape(Nx+1,Ny+1,2)[Nx-2])):.3e}")

# ---- nodal derivative fields (central diff; one-sided at edges) ----
def ddx(f):
    g = np.zeros_like(f)
    g[1:-1,:] = (f[2:,:]-f[:-2,:])/(2*dx)
    g[0,:] = (f[1,:]-f[0,:])/dx
    g[-1,:] = (f[-1,:]-f[-2,:])/dx
    return g
def ddy(f):
    g = np.zeros_like(f)
    g[:,1:-1] = (f[:,2:]-f[:,:-2])/(2*dy)
    g[:,0] = (f[:,1]-f[:,0])/dy
    g[:,-1] = (f[:,-1]-f[:,-2])/dy
    return g

u1x, u1y = ddx(u1), ddy(u1)
u2x, u2y = ddx(u2), ddy(u2)

def interp(field, xp, yp):
    fi = (xp+Lx)/dx; fj = yp/dy
    i = int(np.floor(fi)); j = int(np.floor(fj))
    i = min(max(i,0),Nx-1); j = min(max(j,0),Ny-1)
    tx = fi-i; ty = fj-j
    f00=field[i,j]; f10=field[i+1,j]; f01=field[i,j+1]; f11=field[i+1,j+1]
    return (1-tx)*(1-ty)*f00 + tx*(1-ty)*f10 + (1-tx)*ty*f01 + tx*ty*f11

def J_arc(r0):
    def integrand(th):
        xp = r0*np.cos(th); yp = r0*np.sin(th)
        e11 = interp(u1x,xp,yp); e22 = interp(u2y,xp,yp)
        e12 = 0.5*(interp(u1y,xp,yp)+interp(u2x,xp,yp))
        egg = e11+e22
        w = h*s0*egg + 0.5*h*Ebar*((1-nu)*(e11**2+e22**2+2*e12**2)+nu*egg**2)
        w += 0.5*k*(interp(u1,xp,yp)**2+interp(u2,xp,yp)**2)
        s11 = s0 + Ebar*((1-nu)*e11 + nu*egg)
        s22 = s0 + Ebar*((1-nu)*e22 + nu*egg)
        s12 = Ebar*(1-nu)*e12
        nx = np.cos(th); ny = np.sin(th)
        T1 = h*(s11*nx + s12*ny); T2 = h*(s12*nx + s22*ny)
        du1dx1 = interp(u1x,xp,yp); du2dx1 = interp(u2x,xp,yp)
        return (w*nx - (T1*du1dx1 + T2*du2dx1))*r0
    val,_ = quad(integrand, 0.0, np.pi, limit=400)
    return val

G_analytic = (s0**2*ell/Ebar)*np.tanh(H/(2*ell))
print(f"analytic G_A (H=1,ell=1,sigma0=Ebar=1) = {G_analytic:.8f}")
print(f"analytic gamma_A = {G_analytic:.8f}")
print("r0    J_arc      G_num=2*J_arc   rel.err vs analytic")
for r0 in [0.2,0.3,0.4,0.5,0.6]:
    Ja = J_arc(r0)
    Gnum = 2*Ja/h
    print(f"{r0:4.2f}  {Ja: .6e}   {Gnum: .6e}   {abs(Gnum-G_analytic)/G_analytic:.3e}")
