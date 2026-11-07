import numpy as np
from numpy.polynomial.legendre import leggauss

# ---- mapping t in [-1,1] -> s in [0,inf): s = L*(1+t)/(1-t)
def s_of_t(t,L):
    return L*(1+t)/(1-t)
def dsdt(t,L):
    return 2*L/(1-t)**2

def Kfun(u):
    u=np.asarray(u,dtype=float)
    out=np.zeros_like(u)
    m=np.abs(u)>1e-12
    out[m]=1/np.sinh(u[m]/2)
    return out

def clenshaw_curtis(N):
    # N+1 Lobatto nodes t_j=cos(j pi/N), weights for integral_{-1}^1 f dt
    t=np.cos(np.pi*np.arange(N+1)/N)
    w=np.zeros(N+1)
    for j in range(N+1):
        s=0.0
        for k in range(0,N//2+1):
            b=1.0 if k==0 else 2.0
            if 2*k==N: b=1.0
            s+=b/(4*k*k-1)*np.cos(2*k*j*np.pi/N)
        w[j]=s*2.0/N
    return t,w

def solve_S(eps,a,L=15.0,N=120):
    t,w=clenshaw_curtis(N)
    s=s_of_t(t,L); J=dsdt(t,L)
    # unknown sigma at nodes, with sigma(-1)=1-a, sigma(1)=0
    # collocate at interior nodes i=1..N-1
    # kernel in t: K(t,tau)=J(t)/sinh((s(t)-s(tau))/2)
    # singular part: 2/(t-tau)
    def Kfull(ti,tj):
        return J[j]*Kfun(s[ti]-s[tj]) if False else None
    A=np.zeros((N-1,N-1))
    b=np.zeros(N-1)
    for ii in range(1,N):
        tau=t[ii]
        # RHS = -2 a ln|tanh(X/2)| - int_{-eps}^0 K(X-s) ds, X=s(tau)
        X=s[ii]
        # constant a part: int_0^inf K(X-s) ds = -2 ln|tanh(X/2)|
        const_int=-2*np.log(np.abs(np.tanh(X/2)))
        # known zone integral, Gauss-Legendre on s in [-eps,0]
        xg,wg=leggauss(20)
        ss=-eps/2+eps/2*xg
        zone=- (eps/2)*np.sum(wg*Kfun(X-ss))
        b[ii-1]=-2*a*const_int + zone
        for jj in range(1,N):
            tj=t[jj]
            if ii==jj:
                # 2 * PV int sigma/(t-tau): diagonal after singularity subtraction
                # regular part integral uses (sigma(t)-sigma(tau))/(t-tau); its contribution
                # is approximated by quadrature; the tau-dependent analytic part is
                # 2 sigma(tau) ln|(1-tau)/(1+tau)|
                A[ii-1,jj-1]=0.0  # filled below
            else:
                u=s[ii]-s[jj]
                Kval=J[jj]*Kfun(u)
                # subtract 2/(t_j - t_i)
                Kreg=Kval-2.0/(tj-tau)
                A[ii-1,jj-1]=w[jj]*Kreg
        # diagonal: the regularized integral of (sigma(t)-sigma(tau))/(t-tau) vanishes at t=tau,
        # so the quadrature sum of w_j*Kreg only; plus 2 sigma(tau) ln...
        A[ii-1,ii-1]=2.0*np.log(np.abs((1-tau)/(1+tau)))
        # add endpoint contribution sigma(-1)=1-a to the regular integral?
        # The regular integral approximation uses values at all nodes including endpoints;
        # we already include j=0 in the loop? loop starts at 1, so add endpoint j=0.
        # Add j=0 contribution to all rows:
    # add endpoint j=0 (sigma=1-a) contribution to A and b (it is known, move to RHS)
    rhs=b.copy()
    for ii in range(1,N):
        tau=t[ii]
        u=s[ii]-s[0]
        Kval=J[0]*Kfun(u)
        Kreg=Kval-2.0/(t[0]-tau)
        rhs[ii-1]-=w[0]*Kreg*(1-a)
        # diagonal? no
    return t,w,s,J,A,rhs

def U_value(X,eps,a,t,w,s,J,sigma,L=15.0):
    # U(X) = -int_{-eps}^0 g(X-s) ds - int_0^inf g(X-s) S_+(s) ds
    # g(u)=1/pi^2 ln|coth(u/4)|
    def g(u):
        return (1/np.pi**2)*np.log(np.abs(1/np.tanh(u/4)))
    xg,wg=leggauss(30)
    ss=-eps/2+eps/2*xg
    val=-(eps/2)*np.sum(wg*g(X-ss))
    Splus=a+sigma
    # Splus is defined at t nodes; sigma(-1)=1-a -> Splus(-1)=1, sigma(1)=0 -> Splus(1)=a
    val-=np.sum(w*J*g(X-s)*Splus)
    return val

if __name__=='__main__':
    # test with a from WH a(eps)
    import wh
    for eps in [0.1,0.136,0.2,0.3]:
        a=(-wh.I1(eps)/wh.I2()).real
        t,w,s,J,A,rhs=solve_S(eps,a,L=15,N=120)
        sigma=np.linalg.solve(A,rhs)
        sig=np.concatenate([[1-a],sigma,[0.0]])
        Splus=a+sig
        print('eps',eps,'a',a,'S+(0)',Splus[0],'S+(L)',Splus[-1],'min',Splus.min(),'max',Splus.max())
        print('  U(-eps)',U_value(-eps,eps,a,t,w,s,J,sig))
