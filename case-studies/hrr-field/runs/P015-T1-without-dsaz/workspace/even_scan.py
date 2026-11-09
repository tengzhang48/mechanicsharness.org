import numpy as np, hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import sys
print("even-solution scan (phi'''(0)=0): n, t*, A(pi), C(pi) at A-root")
for nn in [1.01,1.1,1.3,1.5,2.0,3.0,4.0,5.0,6.5]:
    p=1/(nn+1);k=2-p
    def AC(t):
        s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=nn,p=p,k=k),[0,np.pi],[1,0,t,0],
                    rtol=1e-9,atol=1e-11,method='Radau',max_step=0.01)
        phi,d1,d2,d3=s.y[:,-1]
        return k*phi+d2,-(k-1)*d1
    ts=np.linspace(-2.5,1.5,25);vs=[AC(t)[0] for t in ts];roots=[]
    for i in range(len(ts)-1):
        if vs[i]*vs[i+1]<0: roots.append(brentq(lambda t:AC(t)[0],ts[i],ts[i+1],xtol=1e-9))
    if not roots: print("  n=%.2f no A-root"%nn)
    for r in roots:
        A,C=AC(r)
        print("  n=%.2f t*=%.4f A=%.2e C=%.6f"%(nn,r,A,C)); sys.stdout.flush()