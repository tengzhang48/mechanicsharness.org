import numpy as np
from scipy.optimize import brentq
from notch import integrate_many

n=4.0; alpha=2*np.pi/3
ps=np.linspace(-1.6,-0.9,200)

def f_fp(s):
    f,fp=integrate_many(n,s,ps,alpha,nsteps=2000)
    return f,fp

# for each s, find p where f=0 by sign change, then evaluate fp there
def g(s):
    f,fp=f_fp(s)
    for i in range(len(ps)-1):
        if np.isfinite(f[i]) and np.isfinite(f[i+1]) and f[i]*f[i+1]<0:
            # linear interp for p, then local refine
            p0=ps[i]-f[i]*(ps[i+1]-ps[i])/(f[i+1]-f[i])
            f2,fp2=integrate_many(n,s,np.array([p0]),alpha,nsteps=6000)
            return float(fp2[0]), p0
    return np.nan,np.nan

svals=np.linspace(1.835,1.855,9)
vals=[]
for s in svals:
    v,p=g(s); vals.append(v)
    print(f's={s:.5f}  p_f=0 at {p:.6f}  fp(alpha)={v:+.4e}')
for i in range(len(svals)-1):
    if vals[i]*vals[i+1]<0:
        sr=brentq(lambda s:g(s)[0],svals[i],svals[i+1],xtol=1e-10)
        pr=g(sr)[1]
        print('INDEP ROOT: s=%.10f  p_fpp0=%.10f  p=2-s=%.10f'%(sr,pr,2-sr))
