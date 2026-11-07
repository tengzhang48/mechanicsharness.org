import numpy as np
from scipy import special, integrate
from wh import A, Mminus, Mplus, Shat, I1, I2, dA, dMminus

def coefficients(eps, a, N=60):
    an = np.arange(N)+0.5
    qp = 1j*an
    c = np.array([Shat(q,eps,a)/dMminus(q) for q in qp])
    return an, qp, c

def U_full(X, eps, a, V, lam, N=80):
    an, qp, c = coefficients(eps,a,N)
    def Hplus(q):
        return Shat(q,eps,a)/Mminus(q) - np.sum(c/(q-qp))
    qzero = [-1j*(n+0.5) for n in range(N)]
    if V>0:
        qzero.append(-1j/V)
    def dMplus(q):
        h=1e-6
        return (Mplus(q+h,V,lam)-Mplus(q-h,V,lam))/(2*h)
    # residue at q=0 of U=-Hplus/Mplus : Hplus ~ -i a/q, U ~ i a/q => residue i a
    val = a  # residue at 0: U~i a/q, Res=i a, -i*Res = a
    for qj in qzero:
        val += -1j*(-Hplus(qj)/dMplus(qj))*np.exp(1j*qj*X)
    return val

if __name__=='__main__':
    lam=0.1; dbar=0.2
    for V in [0.0]:
        for eps in [0.05,0.1,0.141,0.2,0.3,0.5,1.0]:
            a=-I1(eps)/I2()
            print('V',V,'eps',eps,'a',a,'U',U_full(-eps,eps,a,V,lam),'target',dbar/2)
