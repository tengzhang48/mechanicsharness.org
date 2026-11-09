import numpy as np, time
from notch import integrate_many
ps = np.linspace(-6.0, 2.0, 401)
t0=time.time()
f,fp = integrate_many(4.0, 1.5, ps, 2*np.pi/3, nsteps=3000)
print('one call: %.3fs'%(time.time()-t0))
print('sample f:', f[:3], f[-3:])
