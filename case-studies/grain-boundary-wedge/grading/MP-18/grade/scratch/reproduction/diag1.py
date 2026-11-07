import numpy as np, march

print("N    #zero   top nonzero real parts        log(N)/N   (-re)*N/logN")
for N in [31, 63, 127, 255, 511]:
    M, D1m, D2m, L = march.operator(N)
    ev = np.linalg.eigvals(L)
    re = ev.real
    nz = re[np.abs(re) > 1e-6]
    nzs = np.sort(nz)[::-1]
    nzcount = np.sum(np.abs(re) <= 1e-6)
    c = -nzs[0] * N / np.log(N)
    print("N=%4d nz0=%d  reals=%s  lnN/N=%.5f  c=%.4f" %
          (N, nzcount, np.round(nzs[:3], 5), np.log(N) / N, c))