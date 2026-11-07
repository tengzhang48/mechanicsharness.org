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

def solve_global(a, eps, V, lam=0.1, L=12.0, N1=80, N2=240, N3=120, Ny=150, Xin=None):
    """Monolithic linear solve of moving-frame viscoelastic strip problem.
    Unknowns: U interior, E2 on X<0, T2 on X>=0.
    Constitutive: E2 - V E2_X = T2 - lam V T2_X
    """
    X = make_grid(L, eps, N1, N2, N3) if Xin is None else np.asarray(Xin, float)
    Nx = len(X); dY = np.pi/Ny
    ax = np.zeros(Nx); bx = np.zeros(Nx); cx = np.zeros(Nx)
    for i in range(1, Nx-1):
        hm = X[i]-X[i-1]; hp = X[i+1]-X[i]
        ax[i] = 2.0/((hm+hp)*hm)
        cx[i] = 2.0/((hm+hp)*hp)
        bx[i] = -(ax[i]+cx[i])
    Imin = [i for i in range(1, Nx-1) if X[i] < 0]
    Iplus = [i for i in range(1, Nx-1) if X[i] >= 0]
    imap = {i:k for k,i in enumerate(Imin)}
    pmap = {i:k for k,i in enumerate(Iplus)}
    nu = (Ny-1)*(Nx-2)
    eoff = nu
    toff = eoff + len(Imin)
    N = toff + len(Iplus)
    def ui(i,j): return (j-1)*(Nx-2)+(i-1)
    def ei(i): return eoff+imap[i]
    def ti(i): return toff+pmap[i]
    rows=[]; cols=[]; vals=[]; b=np.zeros(N)
    def add(r,c,v):
        rows.append(r); cols.append(c); vals.append(v)
    # ---------------- Laplace ----------------
    for j in range(1, Ny):
        Yj = j*dY
        for i in range(1, Nx-1):
            r = ui(i,j)
            if j == 1:
                if X[i] >= 0:
                    add(r,r,-2.0/dY**2)
                else:
                    add(r,r,-2.0/(3*dY**2))
                    add(r,ei(i),-2.0/(3*np.pi*dY))
                if j+1 <= Ny-1:
                    add(r,ui(i,j+1),(1.0 if X[i]>=0 else 2.0/3)/dY**2)
            elif j == Ny-1:
                add(r,r,-2.0/dY**2)
                add(r,ui(i,j-1),1.0/dY**2)
                b[r] -= a/dY**2
            else:
                add(r,r,-2.0/dY**2)
                add(r,ui(i,j-1),1.0/dY**2)
                add(r,ui(i,j+1),1.0/dY**2)
            add(r,r,bx[i])
            if i-1 >= 1:
                add(r,ui(i-1,j),ax[i])
            else:
                b[r] -= ax[i]*a
            if i+1 <= Nx-2:
                add(r,ui(i+1,j),cx[i])
            else:
                b[r] -= cx[i]*(a*Yj/np.pi)
    # ---------------- ODE along interface ----------------
    def tfix(x):
        return 1.0 if (-eps <= x < 0.0) else 0.0
    nlap = nu
    for i in range(1, Nx-1):
        row = nlap + (i-1)
        h = X[i+1]-X[i]
        cEi = 1.0 + V/h
        cEi1 = -V/h
        cTi = -1.0 - lam*V/h
        cTi1 = lam*V/h
        def put_E2(i_, coef):
            if i_ <= Nx-2:
                if X[i_] < 0:
                    add(row, ei(i_), coef)
                else:
                    add(row, ui(i_,1), coef*np.pi/dY)
            else:
                b[row] -= coef*a
        def put_T2(i_, coef):
            if i_ <= Nx-2:
                if X[i_] < 0:
                    b[row] -= coef*tfix(X[i_])
                else:
                    add(row, ti(i_), coef)
            else:
                b[row] -= coef*a
        put_E2(i, cEi)
        put_E2(i+1, cEi1)
        put_T2(i, cTi)
        put_T2(i+1, cTi1)
    A = sp.csr_matrix((vals,(rows,cols)), shape=(N,N))
    sol = spla.spsolve(A, b)
    U = np.zeros((Ny-1, Nx-2))
    for j in range(1, Ny):
        for i in range(1, Nx-1):
            U[j-1, i-1] = sol[ui(i,j)]
    E2 = np.zeros(Nx)
    for i in Imin:
        E2[i] = sol[ei(i)]
    for i in Iplus:
        E2[i] = np.pi*U[0, i-1]/dY
    T2 = np.zeros(Nx)
    for i in Imin:
        T2[i] = tfix(X[i])
    for i in Iplus:
        T2[i] = sol[ti(i)]
    U0 = np.zeros(Nx)
    for i in Imin:
        U0[i] = (4*U[0, i-1] - U[1, i-1] - 2*dY*E2[i]/np.pi)/3.0
    U0[Iplus[0]] = 0.0
    return dict(X=X, U=U, E2=E2, T2=T2, U0=U0, dY=dY, sol=sol, a=a, eps=eps, V=V)

if __name__ == '__main__':
    for V in [0.0, 1e-6, 0.1]:
        t0=time.time()
        r = solve_global(np.sqrt(0.2), 0.515, V)
        X=r['X']; U0=r['U0']
        ie = np.argmin(np.abs(X+0.515))
        i0 = np.argmin(np.abs(X))
        print('V',V,'time %.1f'%(time.time()-t0))
        print('  X[i0]',X[i0],'U0(-eps)=',U0[ie],' T2(0)=',r['T2'][i0],
              ' T2+1',r['T2'][i0+1],' E2(-eps)=',r['E2'][ie])
        print('  max|T2 fix dev', np.max(np.abs(r['T2'][X<0]-np.where((X[X<0]>=-0.515),1.,0.))))
