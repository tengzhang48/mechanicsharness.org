import numpy as np, march

# Direct march in u (from march.py docstring logic): d_tau u = D2 s, s = 1 - M w, w = D1 u.
# So d_tau u = D2 (1 - M D1 u) = -D2 M D1 u (constant forcing handled by IC)
# decay operator on u:  A_u = -D2 @ M @ D1
# decay operator on s (via s = -M D1 u): substitute
# Actually let us just march and fit exponential decay at late time and compare eigenvalues.

N = 200
M, D1m, D2m, L = march.operator(N)
print("operator() leading (closest to 0 from below) real eig:", np.max(np.linalg.eigvals(L).real))

# Direct u-march
h = 1.0 / (N + 1)
x = np.arange(1, N + 1) * h
u = np.zeros(N)
s = np.ones(N)
dt = 1e-3
Au = -D2m @ M @ D1m          # d_tau u = D2 s ; d_tau s = -M D1 u  => not closed in u alone
# closed: d_tau u = D2 (1 - M D1 u); the constant part: does 1 map? d_tau u = -D2 M D1 u + D2 1
# At long time s->0. Let's just integrate the full u evolution:
# du/dtau = D2 (1 - M (D1 u))   with s=1- M D1 u
vals = []
for step in range(60000):
    w = D1m @ u
    s = 1.0 - M @ w
    u = u + dt * (D2m @ s)
    if step % 1000 == 0:
        vals.append((step * dt, np.linalg.norm(s)))
vals = np.array(vals)
print("march decay (tau, ||s||):")
for t, nrm in vals:
    print("  %.2f  %.6e" % (t, nrm))
# estimate decay rate from last two
lo = vals[-6:]
rate = -np.polyfit(lo[:, 0], np.log(lo[:, 1]), 1)[0]
print("estimated leading decay rate lambda ~", rate)