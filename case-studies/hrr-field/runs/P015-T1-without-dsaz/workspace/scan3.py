"""Wide scan for ICs at n=6.5 (corrected ODE), report minimum-residual region."""
import numpy as np
from scipy.integrate import solve_ivp
exec(open('scan2.py').read().split('nn=4.0')[0])
nn=6.5
p,k,rhs,endres=make(nn)
best=[]
for a in np.linspace(-2.0,1.0,31):
    for b in np.linspace(-1.0,1.0,31):
        r=endres([a,b],1e-7)
        m=np.max(np.abs(r))
        best.append((m,a,b))
    print("a=%.2f done, best so far: %.5f at (%.3f,%.3f)"%(a,min(best)[0],min(best)[1],min(best)[2]))
best.sort()
print("top:")
for m,a,b in best[:8]:
    print("  m=%.5f  a=%.3f b=%.3f"%(m,a,b))