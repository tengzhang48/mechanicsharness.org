"""
Direct solve_bvp with the CORRECT boundary conditions:
  theta=0: only C(0)=0  => phi'(0)=0   (ONE condition)
  theta=pi: A(pi)=0, C(pi)=0           (TWO conditions)
Total 3 conditions for 4th-order system -> 1 free (scale).
Add normalization phi(0)=1 -> but that over-constrains? No: 3 BCs + the ODE is
4th order => need 4 BCs. The 4th is a normalization. But normalization is a
SCALE choice; it does not remove a solution branch. However solve_bvp needs a
square system: 4 BCs. Use:
  [phi(0)-1, phi'(0), A(pi), C(pi)]
This is exactly hrr2.bc but WITHOUT phi'''(0)=0.
"""
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp

n = 6.5
import hrr2 as H
H.n = n; H.p = 1.0/(n+1.0); H.k = 2.0 - H.p
print("p", H.p, "k", H.k)


def bc(ya, yb):
    k = H.k
    return np.array([
        ya[0] - 1.0,          # phi(0)=1 (scale)
        ya[1],                # phi'(0)=0 (C(0)=0)
        k*yb[0] + yb[2],      # A(pi)=0
        yb[1],                # C(pi)=0
    ])


# initial guess from a rough IVP with phi'''(0)!=0
th = np.linspace(0.0, np.pi, 800)
Y0 = [1.0, 0.0, H.p*(H.p-1.0), 0.5]  # rough
s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-9, atol=1e-11,
              method='Radau', dense_output=True, max_step=0.01)
guess = s.sol(th) if s.success else np.zeros((4, th.size))
guess[0] -= guess[0][0] - 1.0
sol = solve_bvp(H.rhs, bc, th, guess, tol=1e-9, max_nodes=40000, verbose=0)
print("success", sol.success, sol.message)
if sol.success:
    phi = sol.y[0]; d1 = sol.y[1]; d2 = sol.y[2]
    A = H.k*phi + d2
    C = -(H.k-1.0)*d1
    print("A(0)=%.6f C(0)=%.2e A(pi)=%.2e C(pi)=%.2e" % (A[0], C[0], A[-1], C[-1]))
    np.savez("hrr_sol_6p5.npz", th=sol.x, Y=sol.y.T)