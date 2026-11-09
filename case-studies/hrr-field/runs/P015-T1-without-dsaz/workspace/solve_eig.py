"""HRR Mode I: treat p as eigenvalue via solve_bvp.
4 states [phi,phi',phi'',phi'''] + 1 param p.
5 BCs: phi'(0)=0, phi'''(0)=0, phi(0)=1, A(pi)=0, C(pi)=0.
Uses corrected hrr_clean.rhs."""
import numpy as np
import hrr_clean as H
from scipy.integrate import solve_bvp

n = H.N

def rhs_bvp(theta, Y, par):
    p = par[0]
    return H.rhs(theta, Y, n=n, p=p, k=2.0-p)

def bc_bvp(ya, yb, par):
    p = par[0]; k = 2.0-p
    return np.array([
        ya[1],               # phi'(0)=0
        ya[3],               # phi'''(0)=0
        ya[0]-1.0,           # normalization phi(0)=1
        k*yb[0]+yb[2],       # A(pi)=0
        -(k-1.0)*yb[1],      # C(pi)=0  (phi'(pi)=0)
    ])

def solve(p0, guess_amp=0.3):
    th = np.linspace(0.0, np.pi, 300)
    y0 = np.zeros((4, th.size))
    y0[0] = 1.0 + guess_amp*np.cos(th)
    y0[1] = -guess_amp*np.sin(th)
    y0[2] = -guess_amp*np.cos(th)
    sol = solve_bvp(rhs_bvp, bc_bvp, th, y0, p=np.array([p0]),
                    tol=1e-9, max_nodes=20000, verbose=0)
    return sol

if __name__ == "__main__":
    for p0 in [0.5, 0.3, 0.2, 1.0/3]:
        sol = solve(p0)
        print(f"p0={p0:.4f} success={sol.success} {sol.message[:60]} "
              f"p={sol.p} nodes={sol.x.size}")
    # target 1/(n+1)
    print("target p =", 1.0/(n+1))