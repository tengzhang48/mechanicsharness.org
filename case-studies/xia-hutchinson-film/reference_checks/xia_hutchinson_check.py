# Independent check of Xia-Hutchinson Sec 4.2 values from the strip ODE (not agents' code).
import numpy as np
from scipy.integrate import solve_bvp, quad
nu=0.3; Eb=1.0; s0=1.0; ell=1.0; h=1.0; k=h*Eb/ell**2
def U_strip(W):
    # solve h Eb u'' - k u = 0, Eb u'(+-W/2) = -s0 numerically, then integrate energy
    y=np.linspace(-W/2,W/2,401)
    f=lambda y,Y: np.vstack([Y[1], k*Y[0]/(h*Eb)])
    bc=lambda a,b: np.array([Eb*a[1]+s0, Eb*b[1]+s0])
    sol=solve_bvp(f,bc,y,np.zeros((2,y.size)),tol=1e-10,max_nodes=200000)
    u=lambda t: sol.sol(t)[0]; up=lambda t: sol.sol(t)[1]
    pre=quad(lambda t: h*s0*up(t),-W/2,W/2,epsabs=1e-13)[0]
    sto=quad(lambda t: 0.5*h*Eb*up(t)**2+0.5*k*u(t)**2,-W/2,W/2,epsabs=1e-13)[0]
    return pre,sto,sol
rows=[]
for lH in (0.5,1,1.5,2):
    H=ell/lH
    pA,sA,_=U_strip(H); p2,s2,_=U_strip(2*H)
    gA=-(pA+sA)/h*Eb/(s0**2*ell); gB=((p2+s2)-2*(pA+sA))/h*Eb/(s0**2*ell)
    rows.append((lH,gA,gB,np.sqrt((1-nu**2)*gA),np.sqrt((1-nu**2)*gB)))
    print('l/H=%.1f gA=%.6f (tanh %.6f) gB=%.6f (closed %.6f) kIA=%.6f kIB=%.6f  split pre/sto=%.4f'%(lH,gA,np.tanh(H/2),gB,2*np.tanh(H/2)-np.tanh(H),*rows[-1][3:],pA/sA))
gA=lambda H: np.tanh(H/2); gB=lambda H: 2*np.tanh(H/2)-np.tanh(H)
print('H=3l: %.4f %.4f | H=l GB/GA=%.4f sigc ratio=%.4f'%(gA(3),gB(3),gB(1)/gA(1),np.sqrt(gA(1)/gB(1))))
print('strip centre s22: W=1 %.4f W=3 %.4f W=8 %.4f'%tuple(1-1/np.cosh(W/2) for W in (1,3,8)))
ds={0.5:(0.7651,0.5642,0.8344,0.7158),1:(0.4654,0.1642,0.6508,0.3865),1.5:(0.3240,0.0612,0.5430,0.2359),2:(0.2469,0.0280,0.4740,0.1598)}
err=[]
for r in rows:
    for a,b in zip(ds[r[0]],r[1:]): err.append(100*(a/b-1))
print('DeepSeek default errors %%: min %.3f max %.3f'%(min(err),max(err)))
