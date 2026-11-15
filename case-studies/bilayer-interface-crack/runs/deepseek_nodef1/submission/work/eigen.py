import numpy as np

# material properties (plane strain)
nu = 1.0/3.0
E1, E2 = 9.0, 1.0
mu1, mu2 = E1/(2*(1+nu)), E2/(2*(1+nu))
kappa = 3 - 4*nu
Estar1, Estar2 = E1/(1-nu**2), E2/(1-nu**2)
nuprime = nu/(1-nu)

# Dundurs parameters
alpha = (Estar1-Estar2)/(Estar1+Estar2)
beta = 0.5*(mu1*(1-2*nu) - mu2*(1-2*nu))/(mu1*(1-nu)+mu2*(1-nu))
eps = (1/(2*np.pi))*np.log((1-beta)/(1+beta))
print("alpha,beta,eps =", alpha, beta, eps)
print("eps*ln10 deg =", eps*np.log(10)*180/np.pi)

# Airy eigenfunction f(theta) = A e^{i s th} + B e^{-i s th} + C e^{i(s+2)th} + D e^{-i(s+2)th}
def basis(th, s):
    return np.array([np.exp(1j*s*th), np.exp(-1j*s*th),
                     np.exp(1j*(s+2)*th), np.exp(-1j*(s+2)*th)])
def dbasis(th, s):
    return np.array([1j*s*np.exp(1j*s*th), -1j*s*np.exp(-1j*s*th),
                     1j*(s+2)*np.exp(1j*(s+2)*th), -1j*(s+2)*np.exp(-1j*(s+2)*th)])
def ddbasis(th, s):
    return np.array([(1j*s)**2*np.exp(1j*s*th), (-1j*s)**2*np.exp(-1j*s*th),
                     (1j*(s+2))**2*np.exp(1j*(s+2)*th), (-1j*(s+2))**2*np.exp(-1j*(s+2)*th)])

def Pvals(coef, th, s, E, nup):
    f  = coef@basis(th,s)
    fp = coef@dbasis(th,s)
    fpp= coef@ddbasis(th,s)
    Prr = (s+2)*f + fpp
    Ptt = (s+2)*(s+1)*f
    Prt = -(s+1)*fp
    Er  = (Prr - nup*Ptt)/E
    Et  = (Ptt - nup*Prr)/E
    Ert = Prt/(2*0)  # placeholder
    return f,fp,Prr,Ptt,Prt,Er,Et

def UrUth(coef, th, s, E, nup, mu):
    f  = coef@basis(th,s)
    fp = coef@dbasis(th,s)
    fpp= coef@ddbasis(th,s)
    Prr = (s+2)*f + fpp
    Ptt = (s+2)*(s+1)*f
    Prt = -(s+1)*fp
    Ur = (Prr - nup*Ptt)/((s+1)*E)
    # derivative of Ur
    dPrr = (s+2)*fp + coef@np.array([ (1j*s)**3*np.exp(1j*s*th), (-1j*s)**3*np.exp(-1j*s*th),
                     (1j*(s+2))**3*np.exp(1j*(s+2)*th), (-1j*(s+2))**3*np.exp(-1j*(s+2)*th)])
    dPtt = (s+2)*(s+1)*fp
    dUr = (dPrr - nup*dPtt)/((s+1)*E)
    Uth = (Prt/mu - dUr)/s
    return Ur, Uth

def matrix(s, E1=Estar1, E2=Estar2, nup1=nuprime, nup2=nuprime, mu1=mu1, mu2=mu2):
    M = np.zeros((8,8), dtype=complex)
    # material 1: coefficients 0:4; material 2: 4:8
    # f1(pi)=0, f1'(pi)=0
    M[0,0:4] = basis(np.pi,s)
    M[1,0:4] = dbasis(np.pi,s)
    # f2(-pi)=0, f2'(-pi)=0
    M[2,4:8] = basis(-np.pi,s)
    M[3,4:8] = dbasis(-np.pi,s)
    # f1(0)=f2(0)
    M[4,0:4] = basis(0,s); M[4,4:8] = -basis(0,s)
    # f1'(0)=f2'(0)
    M[5,0:4] = dbasis(0,s); M[5,4:8] = -dbasis(0,s)
    # Ur continuity
    Ur1,Ut1 = UrUth(np.eye(4)[0],0,s,E1,nup1,mu1) # dummy
    # compute row operators directly
    def rowUr(E,nup,mu):
        # Ur = (Prr - nup Ptt)/((s+1)E); Prr=(s+2)f+f'', Ptt=(s+2)(s+1)f
        return ((s+2)*basis(0,s)+ddbasis(0,s) - nup*(s+2)*(s+1)*basis(0,s))/((s+1)*E)
    def rowUth(E,nup,mu):
        # Uth = (Prt/mu - dUr)/s
        Prt = -(s+1)*dbasis(0,s)
        # dUr = d/dth[(Prr - nup Ptt)/((s+1)E)]
        dPrr = (s+2)*dbasis(0,s) + np.array([ (1j*s)**3, (-1j*s)**3, (1j*(s+2))**3, (-1j*(s+2))**3])
        dPtt = (s+2)*(s+1)*dbasis(0,s)
        dUr = (dPrr - nup*dPtt)/((s+1)*E)
        return (Prt/mu - dUr)/s
    M[6,0:4] = rowUr(E1,nup1,mu1); M[6,4:8] = -rowUr(E2,nup2,mu2)
    M[7,0:4] = rowUth(E1,nup1,mu1); M[7,4:8] = -rowUth(E2,nup2,mu2)
    return M

for s in [ -0.5+1j*eps, -0.5-1j*eps ]:
    M = matrix(s)
    print("s=",s,"det=",np.linalg.det(M), "svals min=", np.linalg.svd(M, compute_uv=False)[-1])
