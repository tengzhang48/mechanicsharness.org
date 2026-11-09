import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

alpha = 2*np.pi/3

def build(n, Mb, Mf, tol):
    prev = d.run(1.0, alpha, 0.3843, Mb)
    for nn in np.arange(1.25, n+1e-9, 0.25):
        thn = np.linspace(0, alpha, Mb)
        Y = np.zeros((4, Mb))
        for k in range(4): Y[k] = np.interp(thn, prev.x, prev.y[k])
        bcw = lambda ya, yb, p, nn=nn: d.bc(ya, yb, p[0], nn)
        rhsw = lambda x, y, p, nn=nn: d.rhs(x, y, p[0], nn)
        prev = d.solve_bvp(rhsw, bcw, thn, Y, p=np.array([prev.p[0]]), tol=tol, max_nodes=40000)
    # refine fine grid
    thn = np.linspace(0, alpha, Mf)
    Y = np.zeros((4, Mf))
    for k in range(4): Y[k] = np.interp(thn, prev.x, prev.y[k])
    bcw = lambda ya, yb, p: d.bc(ya, yb, p[0], n)
    rhsw = lambda x, y, p: d.rhs(x, y, p[0], n)
    sol = d.solve_bvp(rhsw, bcw, thn, Y, p=np.array([prev.p[0]]), tol=tol, max_nodes=60000)
    s = sol.p[0]
    sr, st, _ = d.stresses(sol.y[0], sol.y[1], sol.y[2], sol.y[3], s, n)
    if st[0] < 0: sr, st = -sr, -st
    return s, sr[0]/st[0], sol.status, np.max(sol.rms_residuals), sol.x.size

for (Mb, Mf, tol) in [(200, 400, 1e-10), (400, 800, 1e-11), (500, 1000, 1e-12), (600, 1500, 1e-12)]:
    s, R, st, rms, nd = build(4.0, Mb, Mf, tol)
    print(f"Mb={Mb} Mf={Mf} tol={tol:g} -> p={s:.8f}  R={R:.8f}  status={st} rms={rms:.1e} nodes={nd}", flush=True)