import numpy as np
from numpy.polynomial.legendre import leggauss

def K(u):
    # 1/sinh(u/2), finite part handled separately
    u=np.asarray(u,dtype=float)
    out=np.zeros_like(u)
    m=np.abs(u)>1e-12
    out[m]=1/np.sinh(u[m]/2)
    return out

def Cfun(X,L):
    # PV integral_0^L ds/sinh((X-s)/2) = 2[ln|tanh((X-L)/2)| - ln|tanh(X/2)|]
    return 2*(np.log(np.abs(np.tanh((X-L)/2))) - np.log(np.abs(np.tanh(X/2))))

def bintegral(X,eps):
    # -int_{-eps}^0 ds/sinh((X-s)/2), X>0, smooth
    xg,wg=leggauss(16)
    mid=-eps/2; half=eps/2
    ss=mid+half*xg
    return -half*np.sum(wg[None,:]*K(X[:,None]-ss[None,:]),axis=1)

def solve(eps,L=15.0,N=600):
    X=np.linspace(0,L,N)
    # regular integral over [0,L] with trapezoid on a fine auxiliary grid
    # A_ij = K(X_i-X_j) for j != i; diagonal coefficient handled by C
    # use the regularized form: A_ij = w_j K(X_i-X_j) for j != i,
    # diagonal: C(X_i) + integral of (K(u)-? ) ... We use pointwise S(X_i),
    # and quadrature weights w_j for the regularized integrand (S(s)-S(X_i))K(X_i-s)
    # Simple approach: use trapezoid weights, and for the term j=i set 0.
    w=np.ones(N); w[0]=0.5; w[-1]=0.5; w*= (L/(N-1))
    A=np.zeros((N,N))
    for i in range(N):
        u=X[i]-X
        A[i,:]=w*K(u)
        A[i,i]=0.0
        A[i,i]=Cfun(X[i],L)
    # Note: the regularized integral uses S(s)-S(X_i); the -S(X_i) part integrated against K
    # with trapezoid excluding j=i is approximate; but the exact PV integral is C, so
    # A S = sum_j w_j K_ij (S_j - S_i) + S_i C = sum_j w_j K_ij S_j + S_i(C - sum_j w_j K_ij)
    for i in range(N):
        A[i,i]=Cfun(X[i],L) - np.sum(w*K(X[i]-X)) + 0.0  # K(0)=0 here
        # Actually K(0) is excluded; the sum is finite. The diagonal should be C - sum_{j!=i} w_j K_ij
        # which is what this is.
    b=bintegral(X,eps)
    S=np.linalg.solve(A,b)
    return X,S,A,b

if __name__=='__main__':
    for N in [200,400,800]:
        for eps in [0.1,0.136,0.2,0.3]:
            X,S,A,b=solve(eps,L=15,N=N)
            print('N',N,'eps',eps,'S0',S[0],'SL',S[-1],'min',S.min(),'max',S.max())
