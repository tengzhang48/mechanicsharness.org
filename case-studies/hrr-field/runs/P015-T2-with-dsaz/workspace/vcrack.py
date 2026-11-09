import numpy as np
from scipy.optimize import fsolve
from notch import integrate_rk4
# crack limit check for several n: should give s=(2n+1)/(n+1)
for n in [3.0,4.0,5.0,9.0]:
    s_hrr=(2*n+1)/(n+1)
    # initial guess p_fpp0 near known ~ -0.66 for n=4; roughly -(2n+1)/... use scan-ish
    guess=[s_hrr,-0.6-(n-4)*0.15]
    try:
        x=fsolve(lambda x: np.array(integrate_rk4(n,x[0],x[1],np.pi,nsteps=6000)),guess,xtol=1e-13)
        print(f'n={n}: s_found={x[0]:.10f}  s_HRR={s_hrr:.10f}  diff={abs(x[0]-s_hrr):.2e}  fpp0={x[1]:.6f}')
    except Exception as e:
        print(f'n={n}: fail {e}')
