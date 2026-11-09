import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

alpha = 2*np.pi/3
n_list = list(np.arange(1.0, 4.0001, 0.25)) + [4.0]
M = 500
prev = None
s_guess = 0.3843
for n in n_list:
    if prev is None:
        sol = d.run(n, alpha, s_guess, M)
    else:
        th_new = np.linspace(0, alpha, M)
        Y = np.zeros((4, M))
        for k in range(4):
            Y[k] = np.interp(th_new, prev.x, prev.y[k])
        a0 = np.array([prev.p[0]])
        bcw = lambda ya, yb, p, n=n: d.bc(ya, yb, p[0], n)
        rhsw = lambda x, y, p, n=n: d.rhs(x, y, p[0], n)
        sol = d.solve_bvp(rhsw, bcw, th_new, Y, p=a0, tol=1e-10, max_nodes=10000)
    s = sol.p[0]
    Srr, Stt, Srt = d.stresses(sol.y[0], sol.y[1], sol.y[2], sol.y[3], s, n)
    R = Srr[0]/Stt[0]
    print(f"n={n:4.2f}  s={s:.8f}  R={R:+.6f}  status={sol.status}  rms={np.max(sol.rms_residuals):.2e}", flush=True)
    prev = sol