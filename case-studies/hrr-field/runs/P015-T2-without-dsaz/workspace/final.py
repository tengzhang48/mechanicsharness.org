import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

alpha = 2*np.pi/3
M = 800

def build(n):
    prev = d.run(1.0, alpha, 0.3843, 300)
    for nn in np.arange(1.25, n+1e-9, 0.25):
        thn = np.linspace(0, alpha, 300)
        Y = np.zeros((4, 300))
        for k in range(4): Y[k] = np.interp(thn, prev.x, prev.y[k])
        bcw = lambda ya, yb, p, nn=nn: d.bc(ya, yb, p[0], nn)
        rhsw = lambda x, y, p, nn=nn: d.rhs(x, y, p[0], nn)
        prev = d.solve_bvp(rhsw, bcw, thn, Y, p=np.array([prev.p[0]]), tol=1e-11, max_nodes=20000)
    return prev

sol = build(4.0)
s = sol.p[0]
# refine on fine grid
thn = np.linspace(0, alpha, M)
Y = np.zeros((4, M))
for k in range(4): Y[k] = np.interp(thn, sol.x, sol.y[k])
bcw = lambda ya, yb, p: d.bc(ya, yb, p[0], 4.0)
rhsw = lambda x, y, p: d.rhs(x, y, p[0], 4.0)
sol = d.solve_bvp(rhsw, bcw, thn, Y, p=np.array([s]), tol=1e-12, max_nodes=40000)
s = sol.p[0]
U, Up, V, P = sol.y
Srr, Stt, Srt = d.stresses(U, Up, V, P, s, 4.0)
# opening amplitude positive: ensure Stt(0)>0
if Stt[0] < 0:
    Srr, Stt, Srt = -Srr, -Stt, -Srt
R = Srr[0]/Stt[0]
print("================================================")
print(f"n=4, alpha=2pi/3 (120 deg opening), plane strain")
print(f"  p (=s)              = {s:.6f}")
print(f"  Srr(0)              = {Srr[0]:.6f}")
print(f"  Stt(0)              = {Stt[0]:.6f}")
print(f"  R = Srr(0)/Stt(0)   = {R:.6f}")
print(f"  status={sol.status} rms={np.max(sol.rms_residuals):.2e} nodes={sol.x.size}")
# finite-energy exponent check: W ~ r^{-s(n+1)}, need s(n+1)<2
print(f"  s(n+1) = {s*5:.4f}  (finite energy requires < 2)")
print("================================================")