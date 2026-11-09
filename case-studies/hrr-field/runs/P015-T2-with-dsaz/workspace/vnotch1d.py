import numpy as np
from scipy.optimize import brentq
from notch import integrate_rk4

n=4.0; alpha=2*np.pi/3
def f_at(s,p): return integrate_rk4(n,s,p,alpha,nsteps=6000)[0]
def fp_at(s,p): return integrate_rk4(n,s,p,alpha,nsteps=6000)[1]

# For a given s, find p with f(alpha)=0 (bracket around -1.2..-1.3)
def p_of_s(s):
    # bracket using the continuation knowledge: p in (-1.5,-1.0)
    a,b=-2.0,-0.5
    fa,fb=f_at(s,a),f_at(s,b)
    # need sign change; scan
    ps=np.linspace(a,b,41); vs=[f_at(s,p) for p in ps]
    for i in range(len(ps)-1):
        if np.isfinite(vs[i]) and np.isfinite(vs[i+1]) and vs[i]*vs[i+1]<0:
            return brentq(lambda p: f_at(s,p), ps[i], ps[i+1], xtol=1e-13)
    return np.nan

# g(s) = f'(alpha) at the p solving f=0
def g(s):
    p=p_of_s(s); return fp_at(s,p)

svals=np.linspace(1.83,1.86,13)
gs=[g(s) for s in svals]
for s,v in zip(svals,gs): print(f's={s:.5f}  fp(alpha)={v:+.6e}')
# find zero
for i in range(len(svals)-1):
    if gs[i]*gs[i+1]<0:
        sroot=brentq(g, svals[i], svals[i+1], xtol=1e-12)
        proot=p_of_s(sroot)
        print('ROOT s=%.12f  p_fpp0=%.12f  p=2-s=%.12f'%(sroot,proot,2-sroot))
