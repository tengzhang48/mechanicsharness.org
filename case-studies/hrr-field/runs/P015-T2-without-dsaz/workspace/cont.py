import numpy as np
from scipy.integrate import solve_bvp
import numpy as np
import bvp as B   # reuse sympy-lambdified RHS/bc builders

# mode-I even linear (n=1) eigenfunction, ODE f''''+(4-2s^2)f''+s^2(2-s)^2 f=0
# f = A cos((2-s)th) + B cos(s th); crack alpha=pi, s=1/2, B=-3A
def lin_seed(s, alpha, th):
    A = 1.0
    Bc = -A*(2 - s)/s        # from A(2-s)+B s = 0
    return A*np.cos((2 - s)*th) + Bc*np.cos(s*th)

def cont(n_list, alpha, s0, seed_th_func, N=300):
    th = np.linspace(0, alpha, N)
    f = seed_th_func(th)
    # derivatives of seed
    h = th[1]-th[0]
    f1 = np.gradient(f, th); f2 = np.gradient(f1, th); f3 = np.gradient(f2, th)
    V = np.zeros_like(f)
    y = np.vstack([f, f1, f2, f3, V])
    p = np.array([s0])
    sol = None
    for n in n_list:
        sol = solve_bvp(B.make_fun(n), B.make_bc(n), th, y, p=p,
                        tol=1e-8, max_nodes=4000, verbose=0)
        th, y, p = sol.x, sol.y, sol.p
        print(f"   n={n}: success={sol.success} s={p[0]:.8f} rms={sol.rms_residuals.max():.2e}", flush=True)
    return sol, th, y, p

if __name__ == "__main__":
    # 1) crack alpha=pi, start near n=1 with s=0.5
    n_list = list(np.linspace(1.0001, 4.0, 25))
    print("crack alpha=pi, continuation in n from 1 to 4 (expect s->0.2):")
    sol, th, y, p = cont(n_list, np.pi, 0.5, lambda t: lin_seed(0.5, np.pi, t))
    print("   crack s =", p[0])
    print("target alpha=2pi/3 continuation in alpha:")
    a_list = list(np.linspace(np.pi, 2*np.pi/3, 30))
    s_cur = p[0]
    seed = lambda t: np.interp(t, th, y[0])
    cur_th, cur_y = th, y
    for a in a_list:
        th = np.linspace(0, a, len(cur_th))
        f = np.interp(th/cur_th[-1]*cur_th[-1], cur_th, cur_y[0])
        f1 = np.gradient(f, th); f2=np.gradient(f1,th); f3=np.gradient(f2,th)
        y0 = np.vstack([f,f1,f2,f3,np.zeros_like(f)])
        sol = solve_bvp(B.make_fun(4.0), B.make_bc(4.0), th, y0, p=np.array([s_cur]),
                        tol=1e-8, max_nodes=4000, verbose=0)
        cur_th, cur_y, s_cur = sol.x, sol.y, sol.p[0]
        print(f"   alpha={a:.5f}: success={sol.success} s={s_cur:.8f}", flush=True)
    print("FINAL s (target) =", s_cur)