import numpy as np
from dissipation_g import fields_T

def gamma_D_psi(r,V,lam=0.1,delta_bar=0.2):
    if V==0: return 0.0
    E1,E2,T1,T2=fields_T(r,V,lam)
    X=r['X']; Y=r['Y']; Xi=X[1:-1]; Yi=Y[1:-1]
    # psi'
    psi=0.5*(E1**2+E2**2)+0.5*lam/(1-lam)*((T1-E1)**2+(T2-E2)**2)
    # X derivatives of psi, E1,E2 on interior grid
    def dX(F):
        out=np.zeros_like(F)
        for k in range(F.shape[1]):
            i=k+1
            if i-1>=1 and i+1<=len(X)-2:
                hm=X[i]-X[i-1]; hp=X[i+1]-X[i]
                out[:,k]=-hp/(hm*(hm+hp))*F[:,k-1]+(hp-hm)/(hm*hp)*F[:,k]+hm/(hp*(hm+hp))*F[:,k+1]
            elif i-1<1:
                out[:,k]=(F[:,k+1]-F[:,k])/(X[i+1]-X[i])
            else:
                out[:,k]=(F[:,k]-F[:,k-1])/(X[i]-X[i-1])
        return out
    dpsi=dX(psi); dE1=dX(E1); dE2=dX(E2)
    # dimensionless integrand for Gamma_D/(tau_c^2/(G_inf tc)) * h/(pi V)? 
    # phi_D = (V/tc) tau_c^2/G_inf * (dpsi - T_i dE_i)
    gdim=dpsi-(T1*dE1+T2*dE2)
    Ix=np.array([np.trapezoid(gdim[j,:],Xi) for j in range(gdim.shape[0])])
    I=np.trapezoid(Ix,Yi)
    # Gamma_D = (1/v) * (h/pi)^2 * tau_c^2/(G_inf tc) * (V) * I *2
    # (1/v)(h/pi)^2*(tau_c^2/(G_inf tc))*V = (pi tc/(Vh))*(h^2/pi^2)*(...)*V = h tc/pi * tau_c^2/(G_inf tc)= h tau_c^2/(pi G_inf)
    # ratio /Gamma0 = [2 h tau_c^2/(pi G_inf) I]/[delta_bar tau_c^2 h/G_inf] = 2 I/(pi delta_bar)
    return 2*I/(np.pi*delta_bar)
