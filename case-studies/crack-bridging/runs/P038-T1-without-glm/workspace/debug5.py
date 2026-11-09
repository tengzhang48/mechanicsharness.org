import numpy as np
from scipy.integrate import quad
import solve as S

Ep=S.Ep; a=S.a; t=S.t
c = S.c_s
SD = (2*t/np.pi)*S.theta(c)
xa = a/c

# spectral b_n (validated: odd ~ exact step coeffs, even ~ 0)
M=1024
k=np.arange(1,M+1); xi=np.cos(k*np.pi/(M+1)); s2=np.sin(k*np.pi/(M+1))**2
g=np.where(np.abs(xi)>xa, SD-t, SD); F=-(4.0*c/Ep)*g
b=np.zeros(M+1); U0=np.ones(M); U1=2*xi
b[1]=(2.0/(M+1))*np.dot(s2*F,U0); b[2]=(2.0/(M+1))*np.dot(s2*F,U1)
for n in range(3,M+1):
    U2=2*xi*U1-U0; b[n]=(2.0/(M+1))*np.dot(s2*F,U2); U0,U1=U1,U2

def P_series(x):
    xi0=x/c; s=0.0; T0=1.0; T1=xi0; s=b[1]*T1
    for n in range(2,M+1):
        T2=2*xi0*T1-T0; s+=b[n]*T2; T0,T1=T1,T2
    return s

def D_J(x):
    sq=np.sqrt(c*c-a*a)
    reg = quad(lambda s:(np.sqrt(c*c-s*s)-sq)/(x-s), -a, a, limit=400)[0]
    J = sq*np.log((x+a)/(x-a)) + reg
    return -4.0/(Ep*np.sqrt(c*c-x*x))*((SD-t)*x + (t/np.pi)*J)

print("  x      P_series      c*sqrt*D_J")
for x in [a*1.0005, 2.7e-3, 3.5e-3, 5e-3, 6.4e-3]:
    print(f"{x*1e3:6.3f}  {P_series(x):+.6e}  {c*np.sqrt(1-(x/c)**2)*D_J(x):+.6e}")

# delta(a) by direct quadrature of -int_a^c D_J dx  (truth check)
val = quad(lambda x: -D_J(x), a, c, points=[a*1.0001], limit=400)[0]
print("delta(a) = -int D_J :", val*1e6, "um  (target 6.8914)")

# delta(a) by quadrature of -int P_series/(sqrt(1-xi^2)) via xi integral
val2 = quad(lambda xi0: -P_series(c*xi0)/np.sqrt(1-xi0*xi0), xa, 1.0, points=[xa,1.0], limit=400)[0]
print("delta(a) = -int P_series/sqrt :", val2*1e6, "um")

# series formula value
s=0.0; V0=1.0; V1=2*xa; s=b[1]*V0
for n in range(2,M+1):
    V2=2*xa*V1-V0; s+=b[n]*V2/n; V0,V1=V1,V2
print("delta(a) series formula:", -np.sqrt(1-xa*xa)*s*1e6, "um")
