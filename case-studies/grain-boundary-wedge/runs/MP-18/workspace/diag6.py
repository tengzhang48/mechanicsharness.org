import numpy as np, march

# Trace the slowest eigenvalue as function of N, identifying the mode by proximity
# to smooth shape (overlap with linear ramp), to separate genuine modes from artifacts.
print("N      lam(slowest,max real)   lam*N/lnN")
for N in [20, 30, 40, 50, 60, 80, 100, 120]:
    M, D1m, D2m, L = march.operator(N)
    ev = np.linalg.eigvals(L)
    j = np.argmax(ev.real)
    lam = -ev.real[j]
    print("%4d   %+.7f   %.4f" % (N, ev.real[j], lam * N / np.log(N)))