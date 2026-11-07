import numpy as np

def fields_T(r, V, lam=0.1):
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
            h0=Y[1]-Y[0]
            E2[j,:]=np.pi*(U[j,:])/h0  # U0=0 for X>=0; for X<0 use E2 from solver? approximate
        elif jj==Ny-1:
            hm=Y[Ny-1]-Y[Ny-2]; hp=Y[Ny]-Y[Ny-1]
            # derivative at last interior using three points U[Ny-2],U[Ny-1],a
            c_m=-hp/(hm*(hm+hp)); c_0=(hp-hm)/(hm*hp); c_p=hm/(hp*(hm+hp))
            E2[j,:]=np.pi*(c_m*U[j-1,:]+c_0*U[j,:]+c_p*a)
        else:
            hm=Y[jj]-Y[jj-1]; hp=Y[jj+1]-Y[jj]
            E2[j,:]=np.pi*(-hp/(hm*(hm+hp))*U[j-1,:]+(hp-hm)/(hm*hp)*U[j,:]+hm/(hp*(hm+hp))*U[j+1,:])
    # fix E2 at j=1 for X<0: use the actual boundary flux E2 stored in r
    # r['E2'] is on the interface X grid; U at j=1 corresponds to all X interior
    # E2=pi*(U1-U0)/h0; U0 = U1 - h0*E2b/pi -> equals E2b. So set directly:
    Xint=X[1:-1]
    E2[0,:]=np.interp(Xint, r['X'], r['E2'])
    # T from ODE backward in X
    T1=np.zeros_like(U); T2=np.zeros_like(U)
    for j in range(U.shape[0]):
        for comp,E,T in ((0,E1,T1),(1,E2,T2)):
            e=E[j,:]; t=np.zeros_like(e)
            e_right=0.0 if comp==0 else a; t_prev=e_right; e_prev=e_right
            for k in range(len(e)-1,-1,-1):
                i=k+1; h=X[i+1]-X[i]
                R=e[k]-V*(e_prev-e[k])/h
                t[k]=(R+lam*V*t_prev/h)/(1.0+lam*V/h)
                t_prev=t[k]; e_prev=e[k]
            T[j,:]=t
    return E1,E2,T1,T2

def gamma_D_ratio(r,V,lam=0.1,delta_bar=0.2):
    if V==0: return 0.0
    E1,E2,T1,T2=fields_T(r,V,lam)
    X=r['X']; Y=r['Y']; Xi=X[1:-1]; Yi=Y[1:-1]
    integ=(T1-E1)**2+(T2-E2)**2
    Ix=np.array([np.trapezoid(integ[j,:],Xi) for j in range(integ.shape[0])])
    I=np.trapezoid(Ix,Yi)
    pref=2/(np.pi*V*(1-lam)*delta_bar)
    return pref*I
