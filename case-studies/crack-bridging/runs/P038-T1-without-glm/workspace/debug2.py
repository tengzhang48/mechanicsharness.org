import numpy as np
from scipy.integrate import quad
import solve as S

Ep=S.Ep; a=S.a; t=S.t
c = S.c_s
SD = (2*t/np.pi)*S.theta(c)

# spectral b_n
def bn_spec(Sr, c, tt, M=1024):
    xa=a/c
    k=np.arange(1,M+1); xi=np.cos(k*np.pi/(M+1)); s2=np.sin(k*np.pi/(M+1))**2
    g=np.where(np.abs(xi)>xa, Sr-tt, Sr); F=-(4.0*c/Ep)*g
    b=np.zeros(M+1); U0=np.ones(M); U1=2*xi
    b[1]=(2.0/(M+1))*np.dot(s2*F,U0); b[2]=(2.0/(M+1))*np.dot(s2*F,U1)
    for n in range(3,M+1):
        U2=2*xi*U1-U0; b[n]=(2.0/(M+1))*np.dot(s2*F,U2); U0,U1=U1,U2
    return b

b = bn_spec(SD, c, t, 1024)

def D_spec(x):  # dislocation density from spectral
    xi = x/c
    # sum b_n T_n(xi)
    s=0.0; T0=1.0; T1=xi
    s=b[1]*T1
    for n in range(2,len(b)):
        T2=2*xi*T1-T0; s+=b[n]*T2; T0,T1=T1,T2
    return s/(c*np.sqrt(1-xi*xi))

def D_exact(x, Sr):
    sq=np.sqrt(c*c-a*a)
    def J(xx):
        reg = quad(lambda s:(np.sqrt(c*c-s*s)-sq)/(xx-s), -a, a, limit=400)[0]
        return sq*np.log((xx+a)/(xx-a))+reg
    return -4.0/(Ep*np.sqrt(c*c-x*x))*((Sr-t)*x + (t/np.pi)*J(x))

for x in [a*1.001, a*1.1, 2.8e-3, 4e-3, 6e-3]:
    print(f"x={x*1e3:.4f}mm  D_spec={D_spec(x):+.6e}  D_J={D_exact(x,SD):+.6e}")

# partial sums of delta series
def delta_series(xi0, Nmax):
    s=0.0; V0=1.0; V1=2*xi0
    s=b[1]*V0
    for n in range(2,Nmax+1):
        V2=2*xi0*V1-V0; s+=b[n]*V2/n; V0,V1=V1,V2
    return -np.sqrt(1-xi0*xi0)*s
for Nmax in [8, 32, 128, 512, 1024]:
    print("delta(a) partial sum N=",Nmax, delta_series(a/c,Nmax)*1e6)
print("target 6.8914")
