import numpy as np, march

print("N     slowest-decay eig (max real part)   -re*N/logN")
prev = None
for N in [15, 25, 35, 45, 65, 95, 135, 195, 275, 395, 555]:
    M, D1m, D2m, L = march.operator(N)
    ev = np.linalg.eigvals(L)
    j = np.argmax(ev.real)
    lam = -ev.real[j]          # decay rate (positive)
    print("%4d   %+.8f      %.4f" % (N, ev.real[j], lam * N / np.log(N)))