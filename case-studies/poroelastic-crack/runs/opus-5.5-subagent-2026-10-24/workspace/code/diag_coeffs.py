import numpy as np
from crack_solver import Params, solve
P = Params(0.2, 0.4, 0.89)
s = 100.0
o = solve(s, 120, P, width=0.5, npts=12)
a = o['L1']['a']*s; b = o['L1']['b']*s
n = 2*np.arange(120)
for i in [0,1,2,4,8,16,32,48,64,80,100,119]:
    print(f"n={n[i]:4d} a_n={a[i]: .3e} (n+1)a_n={(n[i]+1)*a[i]: .3e}  b_n={b[i]: .3e}")
cs = np.cumsum((n+1)*a)/(1-P.nu)
print('partial K sums at N=20,40,60,80,100,120:', cs[[19,39,59,79,99,119]])
