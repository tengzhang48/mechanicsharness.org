import numpy as np
from notch import integrate
# crack check: s=1.8, alpha=pi. find p with f(pi)=0 AND f'(pi)=0
for s in [1.8]:
    for p in [-2.0,-1.5,-1.0,-0.66,-0.5,0.0]:
        sol=integrate(4.0,s,p,np.pi)
        print(f's={s} p={p}: f(pi)={sol.y[0,-1]:+.4f} fp(pi)={sol.y[1,-1]:+.4f}')
