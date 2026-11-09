import numpy as np
from scipy.optimize import fsolve
from notch import integrate_rk4

n=4.0
def F(x, alpha, nsteps=4000):
    s,p = x
    f,fp = integrate_rk4(n,s,p,alpha,nsteps=nsteps)
    return np.array([f, fp])

# solve crack eigenvalue from guess
x = fsolve(lambda x: F(x, np.pi), [1.8, -0.66], xtol=1e-13)
print('crack alpha=pi:', x, ' f,fp=', F(x,np.pi))

# continue in alpha down to 2pi/3
alphas = np.linspace(np.pi, 2*np.pi/3, 30)
xcur = x.copy()
for a in alphas:
    xcur = fsolve(lambda x: F(x, a), xcur, xtol=1e-13)
    print(f'alpha/pi={a/np.pi:.5f}  s={xcur[0]:.10f}  p_fpp0={xcur[1]:.8f}  miss={np.linalg.norm(F(xcur,a)):.2e}')
