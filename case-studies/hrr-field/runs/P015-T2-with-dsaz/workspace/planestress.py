import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

# plane-stress angular ODE (Hutchinson eq. 13), reuse structure from hrr.py
def quad_ps(n,s):
    M=np.zeros((3,3))
    M[0,0]=s*s*(s*s-3*s+3); M[1,1]=3*(1-s)**2; M[2,2]=1.0
    M[0,2]=M[2,0]=s*(3-s)/2.0
    return M
def make_rhs_ps(n,s):
    k=(n-1.0)/2.0; a=n*(s-2.0); b=a+1.0; M=quad_ps(n,s)
    def rhs(x0,x1,x2,x3):
        z=np.array([x0,x1,x2]); v=np.array([x1,x2,x3])
        Q=z@M@z; dQ=2.0*(M@z); dQv=dQ@v; vTMv=v@M@v
        rho=Q**k; g=k*Q**(k-1.0)*dQ; g0,g1,g2=g
        R0=(k*Q**(k-1.0)*(2.0*vTMv)+k*(k-1.0)*Q**(k-2.0)*dQv**2+g0*x2+g1*x3)
        P=g@v
        A=s*(s-3.0)*x0-2.0*x2; Ap=s*(s-3.0)*x1-2.0*x3
        C=s*(2.0*s-3.0)*x0-x2
        coef=2.0*rho-g2*A
        rest=(a*rho*A-R0*A-2.0*P*Ap-rho*s*(s-3.0)*x2+a*b*rho*C+6.0*b*(s-1.0)*(P*x1+rho*x2))
        return -rest/coef
    return rhs

def shoot(n,alpha,guess):
    def F(x):
        sol=solve_ivp(lambda t,y:[y[1],y[2],y[3],make_rhs_ps(n,x[0])(y[0],y[1],y[2],y[3])],
                      [0.0,alpha],[1.0,0.0,x[1],0.0],rtol=1e-12,atol=1e-14,max_step=1e-3)
        return np.array([sol.y[0,-1],sol.y[1,-1]])
    return fsolve(F,guess,xtol=1e-12)

n=4.0
# crack plane stress limit should be s=(2n+1)/(n+1)=1.8
xc=shoot(n,np.pi,[1.8,-0.35])
print('PS crack s=%.8f (expect 1.8)'%xc[0])
# notch
xn=shoot(n,2*np.pi/3,xc)
print('PS notch s=%.8f  p=2-s=%.8f'%(xn[0],2-xn[0]))
