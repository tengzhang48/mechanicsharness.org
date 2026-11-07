import numpy as np
from gbcore import *
from numpy.polynomial.legendre import leggauss

N=40
xg = cheb_lobatto_01(N)
w = bary_weights(xg)
xi, wi = leggauss(8*N); xi=0.5*(xi+1); wi=0.5*wi
L = lagrange_eval(xg, w, xi)
print("partition of unity err:", np.max(np.abs(L.sum(axis=1)-1)))
print("reproduce delta at nodes:", end=" ")
Lex = lagrange_eval(xg, w, xg[1:])
print(np.max(np.abs(Lex - np.eye(N+1)[1:])))
x=0.5
# direct: A[1](x) = ln(x/(1-x)) + int R(x,xi) dxi
numR = np.sum(wi*Rfun(x, xi))
ex = A_exact_const(x)
print("int R direct:", numR, " exact:", ex)
print("ln term:", np.log(x/(1-x)))
