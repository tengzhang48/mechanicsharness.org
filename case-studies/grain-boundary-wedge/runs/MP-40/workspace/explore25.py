exec(open('<round-runs>/runs/MP-40/explore23.py').read().split('for D in')[0])

import numpy as np
D = 24
lam, nodes, Z, On, Vn = run(D)
print("eigen lambda:", ["%.5f" % v for v in lam[:5]])

# rebuild sigma = sum_k c_k T_k(2x-1) from coeffs; Z is in chebyshev basis after run()
# run() returns Z and On,Vn already in chebyshev basis. Need eigenvector:
A = On @ Z; B = Vn @ Z
ev, W = np.linalg.eig(np.linalg.solve(B, A))
cands = [(-ev[k].real, k) for k in range(len(ev)) if abs(ev[k].imag) < 1e-5*max(1, abs(ev[k])) and -ev[k].real > 0]
cands.sort()
k1 = cands[0][1]
a = np.real(W[:, k1])
c = Z @ a


def Tk(k, y):
    return np.cos(k*np.arccos(np.clip(y, -1, 1)))


def sig(x):
    y = 2*np.asarray(x)-1
    return sum(c[k]*Tk(k, y) for k in range(D+1))


f1 = float(sig(np.array([1.0]))[0])
print("lambda1=%.6f  f(0.5)/f(1)=%.6f  f(0.25)/f(1)=%.6f  f(0.75)/f(1)=%.6f"
      % (cands[0][0], float(sig(0.5))/f1, float(sig(0.25))/f1, float(sig(0.75))/f1))

# ---- time stepping: sigma_tau = B^{-1} A a  (backward Euler) ----
L = np.linalg.solve(B, A)
n = L.shape[0]
dt = 0.02
Mbe = np.linalg.inv(np.eye(n) - dt*L)
sigma0 = np.ones(len(nodes))
a = np.linalg.solve(B, sigma0)
tau = 0.0
hist = []
aLate = None
for step in range(200):
    a = Mbe @ a
    tau += dt
    sv = Vn @ (Z @ a)
    hist.append((tau, np.linalg.norm(sv)))
    if aLate is None and tau >= 0.5:
        aLate = a.copy()
hist = np.array(hist)
for lo, hi in [(0.15, 0.4), (0.25, 0.5), (0.3, 0.6), (0.4, 0.8)]:
    mask = (hist[:, 0] > lo) & (hist[:, 0] < hi)
    p = np.polyfit(hist[mask, 0], np.log(hist[mask, 1]), 1)
    print("window [%.2f,%.2f]  decay rate = %.5f" % (lo, hi, -p[0]))
c_late = Z @ aLate
sig_late = lambda x: sum(c_late[k]*Tk(k, 2*np.asarray(x)-1) for k in range(D+1))
print("late profile f(0.5)/f(1)=%.6f  (tau=0.5)" % (float(sig_late(0.5))/float(sig_late(1.0))))