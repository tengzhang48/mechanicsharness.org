import numpy as np, sys
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem, material

p = WavyFrontProblem(A=0.0, harm=1, shift=0.05, Nxf=4096, Nz=64, M=18,
                     Nzmode=2, Lcr=2.0, Lx=3.5, taper=(0.30,0.95))
p.solve()
b = p.total_density()
Ep, mu = p.Ep, p.mu
s2p = np.sqrt(2*np.pi)
C0 = np.array([8*p.K20/(Ep*s2p), 8*p.K10/(Ep*s2p), 4*p.K30/(mu*s2p)])

jr = 0
a0 = 0.05
exact = lambda x: np.where(x < a0, C0*np.sqrt(np.maximum(a0-x,0)), 0.0)
print(" z-row 0: solved vs exact density (y-comp) along x")
for xq in [0.045, 0.03, 0.01, 0.002, -0.002, -0.01, -0.05, -0.15, -0.4, -0.8, -1.2, -1.6, -1.85]:
    i = np.argmin(np.abs(p.xg-xq))
    ex = exact(p.xg[i])
    print(f"x={p.xg[i]:8.4f}  solved {b[1][i,jr]:9.4f}  exact {ex[1]:9.4f} "
          f"ratio {b[1][i,jr]/ex[1] if ex[1]!=0 else np.nan:8.4f}")

# residual pattern of collocation equations
ii, jj = p.collocation()
xc = p.X[ii,jj]
neg = xc<0
rhs = np.zeros((3,len(ii)))
from crack3d import smootherstep
for c in range(3):
    rhs[c] = np.interp(xc, p.xg, p.t_win[c])
    pos = xc>0
    rhs[c][pos] -= (p.K10 if c==1 else (p.K20 if c==0 else p.K30))/np.sqrt(2*np.pi*xc[pos])
# recompute LHS from solved c
H,u = p._window_u()
ns = np.arange(-p.Nzm, p.Nzm+1)
G = np.zeros(p.X.shape, dtype=complex)
for m in range(p.M):
    Tm = np.cos(m*np.arccos(u))
    for i_n,n in enumerate(ns):
        G += p.c[1,m,i_n]*Tm*np.exp(1j*n*p.theta)
gb = [np.zeros_like(G) for _ in range(3)]
gb[1] = G*H
tx,ty,tz = p.op.apply(*gb)
res = np.concatenate([ty[ii,jj].real-rhs[1]])
print(f"\n y-eq residual: rms {np.sqrt((res**2).mean()):.4f}, max {np.abs(res).max():.4f}")
for xr in [(0.0,0.05),(-0.05,0.0),(-0.3,-0.05),(-1.0,-0.3),(-2.0,-1.0)]:
    m = (xc>=xr[0])&(xc<xr[1])
    if m.sum(): print(f"  x in ({xr[0]},{xr[1]}): n={m.sum()}, rms res {np.sqrt((res[m]**2).mean()):.4f}")
