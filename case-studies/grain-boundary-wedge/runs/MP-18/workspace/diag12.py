import numpy as np, march
from march import Smat, D1, D2

# Isolate genuine mode: seed with smooth ramp, project out null space.
# Better: compute eigenpair and check if eigenvector is smooth (physical) or
# a boundary spike (artifact). Look at N=85 zero mode shape.
N = 85
M = Smat(N); D1m = D1(N); D2m = D2(N)
K = -(D2m @ M @ D1m)
ev, V = np.linalg.eig(K)
j0 = np.argmax(ev.real)   # the zero mode
v = V[:, j0].real
x = np.arange(1, N + 1) / (N + 1)
print("zero-mode eig:", ev[j0])
print("shape at x=0.02,0.05,0.1,0.5,1.0:", np.round(np.interp([0.02,0.05,0.1,0.5,1.0],x,v),5))
print("argmax |v| at x=", x[np.argmax(np.abs(v))])

# physical mode: second largest real part
j1 = np.argsort(ev.real)[::-1][1]
v1 = V[:, j1].real
v1 = v1 / v1[-1]
print("\nnext mode eig:", ev[j1])
print("shape( norm f(1)=1) at x=0.1,0.5,1.0:", np.round(np.interp([0.1,0.5,1.0],x,v1),5))