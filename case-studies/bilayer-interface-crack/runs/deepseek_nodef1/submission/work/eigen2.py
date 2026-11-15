import numpy as np
np.set_printoptions(precision=6, suppress=False, linewidth=200)

nu = 1.0/3.0
E1, E2 = 9.0, 1.0
mu1, mu2 = E1/(2*(1+nu)), E2/(2*(1+nu))
Estar1, Estar2 = E1/(1-nu**2), E2/(1-nu**2)
nuprime = nu/(1-nu)
alpha = (Estar1-Estar2)/(Estar1+Estar2)
beta = 0.5*(mu1*(1-2*nu) - mu2*(1-2*nu))/(mu1*(1-nu)+mu2*(1-nu))
eps = (1/(2*np.pi))*np.log((1-beta)/(1+beta))
s = -0.5 + 1j*eps
print("eps=",eps,"s=",s)

def basis(th, s):
    return np.array([np.exp(1j*s*th), np.exp(-1j*s*th),
                     np.exp(1j*(s+2)*th), np.exp(-1j*(s+2)*th)])
def dbasis(th, s):
    return np.array([1j*s*np.exp(1j*s*th), -1j*s*np.exp(-1j*s*th),
                     1j*(s+2)*np.exp(1j*(s+2)*th), -1j*(s+2)*np.exp(-1j*(s+2)*th)])
def ddbasis(th, s):
    return np.array([(1j*s)**2*np.exp(1j*s*th), (-1j*s)**2*np.exp(-1j*s*th),
                     (1j*(s+2))**2*np.exp(1j*(s+2)*th), (-1j*(s+2))**2*np.exp(-1j*(s+2)*th)])
def d3basis(th, s):
    return np.array([(1j*s)**3*np.exp(1j*s*th), (-1j*s)**3*np.exp(-1j*s*th),
                     (1j*(s+2))**3*np.exp(1j*(s+2)*th), (-1j*(s+2))**3*np.exp(-1j*(s+2)*th)])

def matrix(s):
    M = np.zeros((8,8), dtype=complex)
    M[0,0:4] = basis(np.pi,s)
    M[1,0:4] = dbasis(np.pi,s)
    M[2,4:8] = basis(-np.pi,s)
    M[3,4:8] = dbasis(-np.pi,s)
    M[4,0:4] = basis(0,s); M[4,4:8] = -basis(0,s)
    M[5,0:4] = dbasis(0,s); M[5,4:8] = -dbasis(0,s)
    def rowUr(E,nup):
        return ((s+2)*basis(0,s)+ddbasis(0,s) - nup*(s+2)*(s+1)*basis(0,s))/((s+1)*E)
    def rowUth(E,nup,mu):
        Prt = -(s+1)*dbasis(0,s)
        dPrr = (s+2)*dbasis(0,s) + d3basis(0,s)
        dPtt = (s+2)*(s+1)*dbasis(0,s)
        dUr = (dPrr - nup*dPtt)/((s+1)*E)
        return (Prt/mu - dUr)/s
    M[6,0:4] = rowUr(Estar1,nuprime); M[6,4:8] = -rowUr(Estar2,nuprime)
    M[7,0:4] = rowUth(Estar1,nuprime,mu1); M[7,4:8] = -rowUth(Estar2,nuprime,mu2)
    return M

M = matrix(s)
U,S,Vh = np.linalg.svd(M)
coef = Vh[-1].conj()  # nullspace vector
coef = coef/np.linalg.norm(coef)
print("coef=",coef)
print("M@coef=", M@coef)

def field(coef, th, s, E, nup, mu, r=1.0):
    # complex expressions; physical field is real part
    f  = coef@basis(th,s)
    fp = coef@dbasis(th,s)
    fpp= coef@ddbasis(th,s)
    Prr = (s+2)*f + fpp
    Ptt = (s+2)*(s+1)*f
    Prt = -(s+1)*fp
    Ur = (Prr - nup*Ptt)/((s+1)*E)
    dPrr = (s+2)*fp + coef@d3basis(th,s)
    dPtt = (s+2)*(s+1)*fp
    dUr = (dPrr - nup*dPtt)/((s+1)*E)
    Uth = (Prt/mu - dUr)/s
    rs = r**s
    return rs*Prr, rs*Ptt, rs*Prt, r**(s+1)*Ur, r**(s+1)*Uth

# Check T = sig22 + i sig12 at theta=0 for real field
for r in [0.1,1.0,3.0]:
    Prr1,Ptt1,Prt1,ur1,ut1 = field(coef,0,s,Estar1,nuprime,mu1,r)
    Prr2,Ptt2,Prt2,ur2,ut2 = field(coef,0,s,Estar2,nuprime,mu2,r)
    # physical real field: Re of complex. At theta=0, sig22=Re(Ptt), sig12=Re(Prt)
    T = np.real(Ptt1) + 1j*np.real(Prt1)
    print("r",r,"T/r^s =", T/r**s, "Ptt complex", Ptt1, "Prt", Prt1)
    # also physical complex displacement jump? not needed
