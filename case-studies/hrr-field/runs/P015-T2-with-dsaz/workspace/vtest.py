import numpy as np, time
from notch import integrate_many, integrate_rk4, integrate
# validate vectorized vs scalar rk4 at crack
f_sv = integrate_rk4(4.0, 1.8, -0.66, np.pi, nsteps=4000)
f_v, fp_v = integrate_many(4.0, 1.8, np.array([-0.66]), np.pi, nsteps=4000)
print('scalar:', f_sv, ' vec:', f_v[0], fp_v[0])
# speed: full grid s x p
sgrid = np.linspace(1.0, 2.0, 101)
ps = np.linspace(-6.0, 2.0, 401)
t0=time.time()
res = {}
for s in sgrid:
    f, fp = integrate_many(4.0, s, ps, 2*np.pi/3, nsteps=3000)
    res[s]=(f,fp)
print('grid time %.2fs'%(time.time()-t0))
