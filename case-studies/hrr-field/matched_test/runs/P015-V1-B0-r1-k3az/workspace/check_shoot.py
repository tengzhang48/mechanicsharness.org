"""Independent check: shoot the 2nd-order equilibrium ODE directly (no first integral).
h'' = -[nu p^2 h h'^2 + beta q^2 h] / (q^2 + nu h'^2),  q^2 = p^2 h^2 + h'^2,
h(0)=0, h'(0)=1 (scale-invariant ODE). The scaling exponent p=1/(n+1) is correct
iff the trajectory lands on the traction-free condition h'(pi)=0 by itself.
Compare with the first-integral solution from solve.py."""
import numpy as np
from scipy.integrate import solve_ivp
import solve as S

n = S.N
p = 1.0 / (n + 1.0)
nu = (1.0 - n) / n
beta = n / (n + 1.0) ** 2


def ode(t, y):
    h, hp = y
    q2 = p * p * h * h + hp * hp
    hpp = -(nu * p * p * h * hp * hp + beta * q2 * h) / (q2 + nu * hp * hp)
    return [hp, hpp]


sol = solve_ivp(ode, (0.0, np.pi), [0.0, 1.0], method="DOP853",
                dense_output=True, rtol=1e-12, atol=1e-14)
assert sol.success
hpi = sol.sol(np.pi)[1]
print(f"shooting 2nd-order ODE: h'(pi) = {hpi:.3e}  (must be ~0 => exponent consistent)")

tg = np.linspace(0, np.pi, 2001)
h_shoot = sol.sol(tg)[0]
F = S.solve_field(n, S.J_APPLIED)
h_fi = np.interp(tg, F["tg"], F["hg"])
dh = np.max(np.abs(h_shoot - h_fi)) / np.max(np.abs(h_fi))
print(f"max |h_shoot - h_first_integral| / max = {dh:.3e}")

# field values used in the answer, for the record
for tdeg in (0.0, 60.0, 90.0, 120.0, 180.0):
    h, hp, q = S.interp(F, np.radians(tdeg))
    print(f"theta={tdeg:6.1f} deg: h={h:.6f}, h'={hp:.6f}, q={q:.6f}")
assert abs(hpi) < 1e-8 and dh < 1e-8
print("shooting check passed")
