import numpy as np
from solver5 import *

# --- verify J0 derivative closed forms by finite differences ---
x = np.array([0.05, 0.2, 0.4, 0.6, 0.85, 0.999])
J, J1, J2, J3 = J0_pieces(x)
h = 1e-6
JJ = lambda v: J0_pieces(np.array([v]))[0][0]
fd1 = (JJ(x+h)-JJ(x-h))/(2*h)
fd2 = (JJ(x+h)-2*JJ(x)+JJ(x-h))/h**2
fd3 = (JJ(x+2*h)-2*JJ(x+h)+2*JJ(x-h)-JJ(x-2*h))/(2*h**3)
print("J0' err:", np.abs(fd1-J1).max(), " J0'' err:", np.abs(fd2-J2).max(), " J0''' err:", np.abs(fd3-J3).max())

# --- verify Dker derivatives by FD ---
xi = np.linspace(0.05, 0.99, 7)
for k in [1, 2, 3]:
    Dk = Dker(x, xi, k)
    fd = np.zeros_like(Dk)
    hs = 1e-6
    for i, xv in enumerate(x):
        Dp = Dker(np.array([xv+hs]), xi, k-1)[0]; Dm = Dker(np.array([xv-hs]), xi, k-1)[0]
        fd[i] = (Dp-Dm)/(2*hs)
    print(f"D^{k} err:", np.abs(fd-Dk).max())

# --- verify Mm vs quadrature ---
from numpy.polynomial.legendre import leggauss
t, w = leggauss(300); tq=(t+1)/2; wq=w/2
Mm = Mm_polys(8, x)
for m in [2, 5, 8]:
    quad = np.array([np.sum(wq*np.sqrt(1-tq)*(xv-tq)**(m-1)) for xv in x])
    print(f"M_{m} err:", np.abs(Mm[0,:,m-1]-quad).max())

# --- verify qderiv by FD ---
qd = qderiv_vals(x, 6, 4)
for k in [1, 3]:
    qk = lambda v: (qd[k-1, :, :].copy()*0)  # placeholder
for n in [0, 2, 5]:
    qn = lambda v: np.cos(n*np.arccos(2*v-1))
    fd1 = (qn(x+h)-qn(x-h))/(2*h)
    print(f"q_{n}' err:", np.abs(fd1-qd[1,:,n]).max())
