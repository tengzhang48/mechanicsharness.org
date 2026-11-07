import numpy as np

def fields_T_exp(r, V, lam=0.1):
    X=r['X']; Y=r['Y']; U=r['U']; Ny=len(Y)-1; Nx=len(X); a=r['a']
    E1=np.zeros_like(U); E2=np.zeros_like(U)
    for k in range(U.shape[1]):
        i=k+1
        if i-1>=1 and i+1<=Nx-2:
            hm=X[i]-X[i-1]; hp=X[i+1]-X[i]
            E1[:,k]=np.pi*(-hp/(hm*(hm+hp))*U[:,k-1]+(hp-hm)/(hm*hp)*U[:,k]+hm/(hp*(hm+hp))*U[:,k+1])
        elif i-1<1:
            E1[:,k]=np.pi*(U[:,k+1]-U[:,k])/(X[i+1]-X[i])
        else:
            E1[:,k]=np.pi*(U[:,k]-U[:,k-1])/(X[i]-X[i-1])
    for j in range(U.shape[0]):
        jj=j+1
        if jj==1:
            h0=Y[1]-Y[0]; E2[j,:]=np.pi*U[j,:]/h0
        elif jj==Ny-1:
            hm=Y[Ny-1]-Y[Ny-2]; hp=Y[Ny]-Y[Ny-1]
            c_m=-hp/(hm*(hm+hp)); c_0=(hp-hm)/(hm*hp); c_p=hm/(hp*(hm+hp))
            E2[j,:]=np.pi*(c_m*U[j-1,:]+c_0*U[j,:]+c_p*a)
        else:
            hm=Y[jj]-Y[jj-1]; hp=Y[jj+1]-Y[jj]
            E2[j,:]=np.pi*(-hp/(hm*(hm+hp))*U[j-1,:]+(hp-hm)/(hm*hp)*U[j,:]+hm/(hp*(hm+hp))*U[j+1,:])
    Xint=X[1:-1]
    E2[0,:]=np.interp(Xint,r['X'],r['E2'])
    def solve_comp(E,e_right):
        # E shape (nj, n); X full length Nx, E on interior i=1..Nx-2
        # E_X at interior points (nonuniform), one-sided at edges
        EX=np.zeros_like(E)
        for k in range(E.shape[1]):
            i=k+1
            if i-1>=1 and i+1<=Nx-2:
                hm=X[i]-X[i-1]; hp=X[i+1]-X[i]
                EX[:,k]=-hp/(hm*(hm+hp))*E[:,k-1]+(hp-hm)/(hm*hp)*E[:,k]+hm/(hp*(hm+hp))*E[:,k+1]
            elif i-1<1:
                EX[:,k]=(E[:,k+1]-E[:,k])/(X[i+1]-X[i])
            else:
                EX[:,k]=(E[:,k]-E[:,k-1])/(X[i]-X[i-1])
        f=E-V*EX
        T=np.zeros_like(E)
        Tright=e_right*np.ones(E.shape[0])
        for k in range(E.shape[1]-1,-1,-1):
            i=k+1
            h=X[i+1]-X[i]; al=h/(lam*V)
            fr=f[:,k+1] if k+1<f.shape[1] else e_right*np.ones(E.shape[0])
            if al<1e-6:
                loc=0.5*(f[:,k]+fr)*al
            else:
                e=np.exp(-al)
                cf=1.0-(1.0+al)*e
                loc=f[:,k]*(1.0-e)+(fr-f[:,k])*cf/al
            T[:,k]=loc+e*Tright
            Tright=T[:,k]
        return T
    T1=solve_comp(E1,0.0)
    T2=solve_comp(E2,a)
    return E1,E2,T1,T2

def gamma_D_ratio_exp(r,V,lam=0.1,delta_bar=0.2):
    if V==0: return 0.0
    E1,E2,T1,T2=fields_T_exp(r,V,lam)
    X=r['X']; Y=r['Y']; Xi=X[1:-1]; Yi=Y[1:-1]
    integ=(T1-E1)**2+(T2-E2)**2
    Ix=np.array([np.trapezoid(integ[j,:],Xi) for j in range(integ.shape[0])])
    I=np.trapezoid(Ix,Yi)
    pref=2/(np.pi*V*(1-lam)*delta_bar)
    return pref*I
