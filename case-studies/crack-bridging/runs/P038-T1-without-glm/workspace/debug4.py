import numpy as np
import solve as S

Ep=S.Ep; a=S.a; t=S.t
c = S.c_s
SD = (2*t/np.pi)*S.theta(c)
xa = a/c
th_a = np.arccos(xa)

def coeffs(Sr, tt, M=1024):
    k=np.arange(1,M+1); xi=np.cos(k*np.pi/(M+1)); s2=np.sin(k*np.pi/(M+1))**2
    g=np.where(np.abs(xi)>xa, Sr-tt, Sr); F=-(4.0*c/Ep)*g
    b=np.zeros(M+1); U0=np.ones(M); U1=2*xi
    b[1]=(2.0/(M+1))*np.dot(s2*F,U0); b[2]=(2.0/(M+1))*np.dot(s2*F,U1)
    for n in range(3,M+1):
        U2=2*xi*U1-U0; b[n]=(2.0/(M+1))*np.dot(s2*F,U2); U0,U1=U1,U2
    return b

b = coeffs(SD, t)
Fs = -(4*c/Ep)*SD          # F on slit
Fz = -(4*c/Ep)*(SD-t)       # F on zones
dF = Fs-Fz
def c_exact(n):  # exact Chebyshev-U coefficient of the step function F
    if n==1:
        mid = 2*0.5*(xa*np.sqrt(1-xa*xa)+np.arcsin(xa))
        return (2/np.pi)*(dF*mid + Fz*np.pi/2)
    br = lambda th: np.sin((n-1)*th)/(n-1) - np.sin((n+1)*th)/(n+1)
    mid = br(np.pi/2) - br(th_a)
    return (2/np.pi)*dF*mid

print(" n   b_discrete      c_exact(step)")
for n in [1,2,3,4,5,6,7,10,15,20,30,50,100,200]:
    print(f"{n:4d}  {b[n]:+.6e}  {c_exact(n):+.6e}")

# series with EXACT step coefficients
def delta_series_bn(bn, xi0, Nmax):
    s=0.0; V0=1.0; V1=2*xi0
    s=bn[1]*V0
    for n in range(2,Nmax+1):
        V2=2*xi0*V1-V0; s+=bn[n]*V2/n; V0,V1=V1,V2
    return -np.sqrt(1-xi0*xi0)*s
bex = np.zeros(2001)
for n in range(1,2001):
    bex[n] = c_exact(n) if n>1 else c_exact(1)
for Nmax in [8,32,128,512,2000]:
    print("delta(a) with EXACT step coeffs, N=",Nmax, delta_series_bn(bex, xa, Nmax)*1e6)
print("literature Dugdale:", (8*t*a/(np.pi*Ep))*np.log(1/np.cos(th_a))*1e6)
