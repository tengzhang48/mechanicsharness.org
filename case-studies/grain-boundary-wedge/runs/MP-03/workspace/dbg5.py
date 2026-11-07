import numpy as np
from gbcore import cheb_lobatto_01, diff_mat
N=20
xg = cheb_lobatto_01(N); D = diff_mat(xg); D2 = D@D; D3 = D2@D
# test on polynomial p = x^5 - 3x^2 + 2x
p = xg**5 - 3*xg**2 + 2*xg
print("p'  err:", np.max(np.abs(D@p - (5*xg**4 - 6*xg + 2))))
print("p'' err:", np.max(np.abs(D2@p - (20*xg**3 - 6))))
print("p''' err:", np.max(np.abs(D3@p - 60*xg**2)))
