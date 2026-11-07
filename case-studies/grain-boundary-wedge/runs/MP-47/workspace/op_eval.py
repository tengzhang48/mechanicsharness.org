import numpy as np
from explore_eig import clencurt, chebD, Rreg
from op2 import cheb_interp_matrix

def build_Op_eval(x, N, xe, xq, wq):
    """p at Chebyshev nodes -> int S(xe,xi)p(xi)dxi at evaluation nodes xe."""
    Icall=cheb_interp_matrix(x,xq,N)     # (nq, N+1)
    ne=len(xe)
    Op=np.zeros((ne,N+1))
    for ii in range(ne):
        xi=xe[ii]
        coef1=wq/(xi-xq)
        Op[ii,:]+= coef1@Icall
        s=np.sum(coef1)
        # ln term: needs p(xi) which is interp of p at xi. Use interpolation row.
        Irow=cheb_interp_matrix(x,np.array([xi]),N)[0]  # (N+1,)
        Op[ii,:]+= (np.log(abs(xi/(1.0-xi)))-s)*Irow
        coef3=wq*Rreg(xi,xq)
        Op[ii,:]+= coef3@Icall
    return Op