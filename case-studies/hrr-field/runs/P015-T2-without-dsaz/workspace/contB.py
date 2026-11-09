import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location("d", "disp.py"); d = importlib.util.module_from_spec(spec); spec.loader.exec_module(d)

alpha = 2*np.pi/3; M = 500
# start from Mode B at n=4, go DOWN in n
sol = d.run(4.0, alpha, 0.35, M)
nl = list(np.arange(4.0, 0.999, -0.25))
for n in nl[1:]:
    thn = np.linspace(0, alpha, M)
    Y = np.zeros((4, M))
    for k in range(4): Y[k] = np.interp(thn, sol.x, sol.y[k])
    bcw = lambda ya, yb, p, n=n: d.bc(ya, yb, p[0], n)
    rhsw = lambda x, y, p, n=n: d.rhs(x, y, p[0], n)
    sol = d.solve_bvp(rhsw, bcw, thn, Y, p=np.array([sol.p[0]]), tol=1e-10, max_nodes=10000)
    s = sol.p[0]
    Srr, Stt, Srt = d.stresses(sol.y[0], sol.y[1], sol.y[2], sol.y[3], s, n)
    print(f"n={n:4.2f} s={s:.6f} R={Srr[0]/Stt[0]:+.5f} status={sol.status} rms={np.max(sol.rms_residuals):.2e}", flush=True)