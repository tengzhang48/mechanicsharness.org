import numpy as np, march
from march import Smat, D3

# Correct reduced operator L = -M @ D3  (Smat already carries 1/pi).
# Check spectrum across N: is there a genuine discrete leading eigenvalue?
print("N     slowest-decay eig (-re smallest magnitude, nonzero)   top 4 reals")
for N in [15, 25, 35, 45, 65, 95, 135]:
    M = Smat(N)
    D3m = D3(N)
    L = -(M @ D3m)
    ev = np.linalg.eigvals(L)
    re = ev.real
    nz = re[np.abs(re) > 1e-8]
    srt = np.sort(nz)[::-1]  # closest to 0 from below
    print("%4d  %s" % (N, np.round(srt[:4], 6)))