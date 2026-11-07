import numpy as np, time
from numpy.polynomial.legendre import leggauss
t0=time.time()
from cheb2 import Rval, gA, gB, brute_R
print("import",time.time()-t0)
for fam,k in [('A',1),('B',2),('A',3),('B',0)]:
    for Z in [0.2,0.5,0.7]:
        print(fam,k,round(Z,2), "Rval=%.6f"%Rval(k,fam,Z), "brute=%.6f"%brute_R(lambda W:gA(k,W) if fam=='A' else gB(k,W),Z))
print("done",time.time()-t0)