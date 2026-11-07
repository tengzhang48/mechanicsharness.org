import numpy as np, march

N = 60
M, D1m, D2m, L = march.operator(N)
ev, V = np.linalg.eig(L)
j = np.argmin(np.abs(ev))
v = V[:, j].real
v = v / v[-1]
x = np.arange(1, N + 1) / (N + 1)
print("near-zero eig:", ev[j])
print("v sample at x=0.1,0.25,0.5,0.75,0.95:", np.round(np.interp([0.1,0.25,0.5,0.75,0.95], x, v), 5))
print("D2 v norm:", np.linalg.norm(D2m @ v), " D1 D2 v norm:", np.linalg.norm(D1m @ D2m @ v))

# which row of M@D1@D2 breaks? check residuals row by row
r = L @ v
print("residual norm", np.linalg.norm(r))
print("residual first 5:", np.round(r[:5], 8), " last 5:", np.round(r[-5:], 8))

# is 0 exactly an eigenvalue because M is singular?
print("rank/cond of M:", np.linalg.matrix_rank(M), np.linalg.cond(M))
# Check nullspace of M
sM = np.linalg.svd(M)[1]
print("smallest singular values of M:", np.round(sM[-5:], 8))