"""Isolate the optional amplitude claim and an invalid K-validation claim."""
import numpy as np
from scipy.linalg import eig
import solve
import sblib as sb
from gblib import T_moment

coefficients = []
for n in [128, 256, 512]:
    L, x, _ = solve.build_FD_cellcentered(n)
    w, left, right = eig(L, left=True, right=True)
    i = np.argmax(w.real)
    v = right[:,i].real
    v /= v[-1]
    ell = left[:,i].real
    amplitude = (ell @ np.ones(n))/(ell @ v)
    coefficients.append(amplitude)
    rate = -w[i].real
    euler_rate = np.log1p(rate*2e-4)/2e-4
    print('N=',n,'modal_c1=',amplitude,'lambda1=',rate,'Euler_rate=',euler_rate)
print('First-order amplitude extrapolation from N=256,512:',2*coefficients[-1]-coefficients[-2])

wrong = np.array([-0.1])
print('Incorrect density: c0=-0.1')
for x in [.1,.5,.9]:
    s = np.sqrt(1-x)
    traction_integral = 2*wrong[0]*sb.R_all(s,0)[0]
    print('x=',x,'T(vhat)=',traction_integral,'target=1')
print('Wrong-density K_tip=',-np.pi*np.sqrt(2)*wrong[0])
for e in [1e-5,1e-7,1e-9]:
    outside = np.sqrt(2*e)*(1-sb.T_outside(wrong,1+e))
    print('eps=',e,'K_out=',outside)

print('Exact mouth-limit evidence: R_1(s) = T_0(1-s*s)/2.')
for x in [1e-2,1e-4,1e-6,1e-8]:
    print('x=',x,'T0/2=',T_moment(x,0)/2)
