import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def make_grid(L,eps,N1=80,N2=240,N3=120):
    t=np.linspace(0,1,N1+1); u=(1-np.cos(np.pi*t))/2
    x1=-L+(L-eps)*u
    t=np.linspace(0,1,N2+1); u=(1-np.cos(np.pi*t))/2
    x2=-eps+eps*u
    t=np.linspace(0,1,N3+1); u=(1-np.cos(np.pi*t))/2
    x3=L*u
    return np.concatenate([x1,x2[1:],x3[1:]])

def solve_global(a, eps, V, lam=0.1, L=12.0, N1=80, N2=240, N3=120, Ny=150,
                 Xin=None, Yclust=True):
    X = make_grid(L, eps, N1, N2, N3) if Xin is None else np.asarray(Xin, float)
    if Yclust:
        t=np.linspace(0,1,Ny+1); Y=np.pi*(1-np.cos(np.pi*t))/2
    else:
        Y=np.linspace(0,np.pi,Ny+1)
    Nx=len(X)
    ax=np.zeros(Nx); bx=np.zeros(Nx); cx=np.zeros(Nx)
    for i in range(1,Nx-1):
        hm=X[i]-X[i-1]; hp=X[i+1]-X[i]
        ax[i]=2.0/((hm+hp)*hm); cx[i]=2.0/((hm+hp)*hp); bx[i]=-(ax[i]+cx[i])
    Imin=[i for i in range(1,Nx-1) if X[i]<0]
    Iplus=[i for i in range(1,Nx-1) if X[i]>=0]
    imap={i:k for k,i in enumerate(Imin)}; pmap={i:k for k,i in enumerate(Iplus)}
    nu=(Ny-1)*(Nx-2); eoff=nu; toff=eoff+len(Imin); Ntot=toff+len(Iplus)
    def ui(i,j): return (j-1)*(Nx-2)+(i-1)
    def ei(i): return eoff+imap[i]
    def ti(i): return toff+pmap[i]
    rows=[];cols=[];vals=[]; b=np.zeros(Ntot)
    def add(r,c,v): rows.append(r);cols.append(c);vals.append(v)
    def tfix(x): return 1.0 if (-eps<=x<0.0) else 0.0
    # Laplace
    for j in range(1,Ny):
        if j==1:
            hm=Y[1]-Y[0]; hp=Y[2]-Y[1]
        elif j==Ny-1:
            hm=Y[Ny-1]-Y[Ny-2]; hp=Y[Ny]-Y[Ny-1]
        else:
            hm=Y[j]-Y[j-1]; hp=Y[j+1]-Y[j]
        for i in range(1,Nx-1):
            r=ui(i,j)
            # Y part
            if j==1:
                h0=Y[1]-Y[0]; h1=Y[2]-Y[1]
                if X[i]>=0:
                    add(r,r,-2.0/(h0*h1))
                    add(r,ui(i,j+1),2.0/(h1*(h0+h1)))
                else:
                    add(r,r,-2.0/(h1*(h0+h1)))
                    add(r,ui(i,j+1),2.0/(h1*(h0+h1)))
                    add(r,ei(i),-2.0/(np.pi*(h0+h1)))
            elif j==Ny-1:
                c_m=2.0/(hm*(hm+hp)); c_0=-2.0/(hm*hp); c_p=2.0/(hp*(hm+hp))
                add(r,r,c_0)
                add(r,ui(i,j-1),c_m)
                b[r]-=c_p*a
            else:
                c_m=2.0/((hm+hp)*hm); c_0=-(c_m+2.0/((hm+hp)*hp)); c_p=2.0/((hm+hp)*hp)
                add(r,r,c_0)
                add(r,ui(i,j-1),c_m)
                add(r,ui(i,j+1),c_p)
            # X part
            add(r,r,bx[i])
            if i-1>=1: add(r,ui(i-1,j),ax[i])
            else: b[r]-=ax[i]*a
            if i+1<=Nx-2: add(r,ui(i+1,j),cx[i])
            else: b[r]-=cx[i]*(a*Y[j]/np.pi)
    # ODE
    nlap=nu
    for i in range(1,Nx-1):
        row=nlap+(i-1); h=X[i+1]-X[i]
        cEi=1.0+V/h; cEi1=-V/h; cTi=-1.0-lam*V/h; cTi1=lam*V/h
        def put_E2(i_,coef):
            if i_<=Nx-2:
                if X[i_]<0: add(row,ei(i_),coef)
                else: add(row,ui(i_,1),coef*np.pi/(Y[1]-Y[0]))
            else: b[row]-=coef*a
        def put_T2(i_,coef):
            if i_<=Nx-2:
                if X[i_]<0: b[row]-=coef*tfix(X[i_])
                else: add(row,ti(i_),coef)
            else: b[row]-=coef*a
        put_E2(i,cEi); put_E2(i+1,cEi1); put_T2(i,cTi); put_T2(i+1,cTi1)
    A=sp.csr_matrix((vals,(rows,cols)),shape=(Ntot,Ntot))
    sol=spla.spsolve(A,b)
    U=np.zeros((Ny-1,Nx-2))
    for j in range(1,Ny):
        for i in range(1,Nx-1): U[j-1,i-1]=sol[ui(i,j)]
    E2=np.zeros(Nx)
    for i in Imin: E2[i]=sol[ei(i)]
    h0=Y[1]-Y[0]
    for i in Iplus: E2[i]=np.pi*U[0,i-1]/h0
    T2=np.zeros(Nx)
    for i in Imin: T2[i]=tfix(X[i])
    for i in Iplus: T2[i]=sol[ti(i)]
    U0=np.zeros(Nx)
    for i in Imin: U0[i]=U[0,i-1]-h0*E2[i]/np.pi
    U0[Iplus[0]]=0.0
    return dict(X=X,Y=Y,U=U,E2=E2,T2=T2,U0=U0,dY=h0,sol=sol,a=a,eps=eps,V=V,Ny=Ny)
