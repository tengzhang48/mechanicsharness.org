import numpy as np, march
np.set_printoptions(precision=4, suppress=True, linewidth=150)

N = 40
M, D1m, D2m, L = march.operator(N)
ev, V = np.linalg.eig(L)
re = ev.real
# near-zero mode
i0 = np.argmin(np.abs(ev))
print("near-zero eig:", ev[i0])
v0 = V[:, i0].real
print("near-zero eigenvector (first/last 8):")
print(np.round(v0[:8], 5))
print(np.round(v0[-8:], 5))

# physical leading mode (largest |real| but negative - that's fastest decay). 
# We want slowest decay = max real part. Among nonzero:
mask = np.abs(ev.real) > 1e-6
j = np.argmax(np.where(mask, ev.real, -1e9))
print("\nnext eig (max real, nonzero):", ev[j])
vj = V[:, j].real
print("eigvec:", np.round(vj[:8], 5), "...", np.round(vj[-4:], 5))

# Check: is the near-zero eigenvector just a highly oscillatory / boundary artifact?
print("\nnorm v0", np.linalg.norm(v0), " argmax|x*v0|", np.argmax(np.abs(v0)))
x = np.arange(1, N + 1) / (N + 1)
print("v0 at x=0.5 index:", v0[N // 2 - 1])