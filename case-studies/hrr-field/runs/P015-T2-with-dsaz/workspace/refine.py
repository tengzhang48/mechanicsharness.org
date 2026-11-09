import numpy as np
from scipy.optimize import fsolve
from notch import integrate_rk4

n=4.0; alpha=2*np.pi/3
# start from continuation result
x0=np.array([1.8438247397,-1.26307900])
for ns in [1000,2000,4000,8000]:
    x=fsolve(lambda x: np.array(integrate_rk4(n,x[0],x[1],alpha,nsteps=ns)), x0, xtol=1e-14)
    s,pf=x
    print(f'nsteps={ns:5d}  s={s:.12f}  fpp0={pf:.12f}  p=2-s={2-s:.12f}')
