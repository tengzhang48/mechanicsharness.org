import numpy as np
from eig import eigen_s
for (N,M) in [(48,None),(48,2048),(60,None),(60,4096),(72,None),(72,6144),(96,None)]:
    xg, lam, V, lams = eigen_s(N=N, M=M, verbose=False)
    print(f"N={N}, M={M if M else 8*N}: lam1={lams[0]:.6f} lam2={lams[1]:.4f} lam3={lams[2]:.3f}")
