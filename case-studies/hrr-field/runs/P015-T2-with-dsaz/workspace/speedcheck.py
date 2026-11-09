import numpy as np, time
from notch import integrate, integrate_rk4
# crack: s=1.8, p=-0.66, alpha=pi
f_sv = integrate(4.0, 1.8, -0.66, np.pi).y[0,-1]
fp_sv = integrate(4.0, 1.8, -0.66, np.pi).y[1,-1]
t0=time.time(); f_rk, fp_rk = integrate_rk4(4.0, 1.8, -0.66, np.pi, nsteps=4000); t1=time.time()
print('solve_ivp:', f_sv, fp_sv)
print('rk4     :', f_rk, fp_rk, ' dt=%.4fs'%(t1-t0))
# notch alpha
t0=time.time(); fa, fpa = integrate_rk4(4.0, 1.5, -1.0, 2*np.pi/3, nsteps=4000); t1=time.time()
print('notch rk4:', fa, fpa, ' dt=%.4fs'%(t1-t0))
