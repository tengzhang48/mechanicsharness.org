import numpy as np

# ---------------- material constants ----------------
nu = 1.0/3.0
E1, E2 = 9.0, 1.0
mu1, mu2 = E1/(2*(1+nu)), E2/(2*(1+nu))
kappa = 3-4*nu
Estar1, Estar2 = E1/(1-nu**2), E2/(1-nu**2)
nuprime = nu/(1-nu)
alpha = (Estar1-Estar2)/(Estar1+Estar2)
beta = 0.5*(mu1*(1-2*nu) - mu2*(1-2*nu))/(mu1*(1-nu)+mu2*(1-nu))
eps = (1/(2*np.pi))*np.log((1-beta)/(1+beta))
s = -0.5 + 1j*eps

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
    M[0,0:4] = basis(np.pi,s); M[1,0:4] = dbasis(np.pi,s)
    M[2,4:8] = basis(-np.pi,s); M[3,4:8] = dbasis(-np.pi,s)
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

def nullspace(s):
    M = matrix(s)
    U,S,Vh = np.linalg.svd(M)
    c = Vh[-1].conj(); c = c/np.linalg.norm(c)
    return c[:4], c[4:]

c1, c2 = nullspace(s)

# normalized auxiliary fields: A = Re(F), B = Im(F)
# K_A, K_B from theta=0 stress
def polar_complex(c, th, r, E, nup, mu):
    f  = c@basis(th,s); fp = c@dbasis(th,s); fpp= c@ddbasis(th,s)
    Prr = (s+2)*f + fpp
    Ptt = (s+2)*(s+1)*f
    Prt = -(s+1)*fp
    Ur = (Prr - nup*Ptt)/((s+1)*E)
    dPrr = (s+2)*fp + c@d3basis(th,s)
    dPtt = (s+2)*(s+1)*fp
    dUr = (dPrr - nup*dPtt)/((s+1)*E)
    Uth = (Prt/mu - dUr)/s
    rs = r**s
    return rs*Prr, rs*Ptt, rs*Prt, r**(s+1)*Ur, r**(s+1)*Uth

Prr,Ptt,Prt,ur,ut = polar_complex(c1,0,1.0,Estar1,nuprime,mu1)
KA = Ptt                      # complex expression; real field T = Re?  we found T_A = Ptt
KB = 1j*np.conj(Ptt)          # real field B
print("KA =",KA," KB =",KB)
print("check B: Im(Ptt)+i Im(Prt) =", np.imag(Ptt)+1j*np.imag(Prt))

# constants for G = Cc |K|^2
Cc = ((kappa+1)/mu1 + (kappa+1)/mu2)/(16*np.cosh(np.pi*eps)**2)
print("C =",Cc, " |KA|^2*C =", Cc*abs(KA)**2)

# ---------------- analytic Cartesian fields ----------------
def field_cart(mat, c, r, th, derivs=False):
    """complex expressions; physical field is Re or Im.
       returns sig_xx,sig_yy,sig_xy, ux,uy and optionally derivatives dux/dx etc."""
    E, nup, mu = (Estar1,nuprime,mu1) if mat==1 else (Estar2,nuprime,mu2)
    Prr,Ptt,Prt,Ur_ph,Ut_ph = polar_complex(c, th, r, E, nup, mu)
    # Ur_ph = r^{s+1} U_r, Ut_ph = r^{s+1} U_theta
    U_r = Ur_ph/r      # actually U_r(θ) * r^s?  careful
    # Let's compute directly from angular functions:
    f  = c@basis(th,s); fp = c@dbasis(th,s); fpp= c@ddbasis(th,s)
    P_rr = (s+2)*f + fpp; P_tt=(s+2)*(s+1)*f; P_rt=-(s+1)*fp
    Ur = (P_rr - nup*P_tt)/((s+1)*E)
    dPrr = (s+2)*fp + c@d3basis(th,s); dPtt = (s+2)*(s+1)*fp
    dUr = (dPrr - nup*dPtt)/((s+1)*E)
    Uth = (P_rt/mu - dUr)/s
    dUr_dth = dUr
    # Uth derivative: Uth = (Prt/mu - dUr)/s
    dPrt = -(s+1)*ddbasis(th,s)@c
    d2Prr = (s+2)*ddbasis(th,s)@c + np.array([(1j*s)**4,(-1j*s)**4,(1j*(s+2))**4,(-1j*(s+2))**4])@c
    d2Ptt = (s+2)*(s+1)*ddbasis(th,s)@c
    d2Ur = (d2Prr - nup*d2Ptt)/((s+1)*E)
    dUth = (dPrt/mu - d2Ur)/s
    rs = r**s
    Sxx = rs*(P_rr*np.cos(th)**2 + P_tt*np.sin(th)**2 - 2*P_rt*np.sin(th)*np.cos(th))
    Syy = rs*(P_rr*np.sin(th)**2 + P_tt*np.cos(th)**2 + 2*P_rt*np.sin(th)*np.cos(th))
    Sxy = rs*((P_rr-P_tt)*np.sin(th)*np.cos(th) + P_rt*(np.cos(th)**2-np.sin(th)**2))
    rsp1 = r**(s+1)
    ux = rsp1*(Ur*np.cos(th) - Uth*np.sin(th))
    uy = rsp1*(Ur*np.sin(th) + Uth*np.cos(th))
    if not derivs:
        return Sxx,Syy,Sxy,ux,uy
    # derivatives wrt x
    # d/dx = cos th d/dr - sin th /r d/dth
    # ur = r^{s+1}Ur, uth=r^{s+1}Uth
    # build derivatives of ux,uy wrt r and th
    ux_r = (s+1)*r**s*(Ur*np.cos(th)-Uth*np.sin(th))
    uy_r = (s+1)*r**s*(Ur*np.sin(th)+Uth*np.cos(th))
    ux_th = rsp1*(dUr*np.cos(th)-Ur*np.sin(th)-dUth*np.sin(th)-Uth*np.cos(th))
    uy_th = rsp1*(dUr*np.sin(th)+Ur*np.cos(th)+dUth*np.cos(th)-Uth*np.sin(th))
    ux_x = np.cos(th)*ux_r - np.sin(th)/r*ux_th
    uy_x = np.cos(th)*uy_r - np.sin(th)/r*uy_th
    ux_y = np.sin(th)*ux_r + np.cos(th)/r*ux_th
    uy_y = np.sin(th)*uy_r + np.cos(th)/r*uy_th
    return Sxx,Syy,Sxy,ux,uy,ux_x,uy_x,ux_y,uy_y

# check strain compatibility for field A: eps_xx = dux/dx etc vs compliance
def check_disp(mat, r=1.0, th=0.7):
    Sxx,Syy,Sxy,ux,uy,uxx,uyx,uxy,uyy = field_cart(mat,c1,r,th,True)
    # physical A = Re
    exx = np.real(uxx); eyy=np.real(uyy); exy=0.5*np.real(uxy+uyx)
    E,nup,mu = (Estar1,nuprime,mu1) if mat==1 else (Estar2,nuprime,mu2)
    # plane strain compliance
    sxx=np.real(Sxx); syy=np.real(Syy); sxy=np.real(Sxy)
    exx_c = (sxx - nup*syy)/E
    eyy_c = (syy - nup*sxx)/E
    exy_c = sxy/(2*mu)
    return (exx,exx_c,eyy,eyy_c,exy,exy_c)

print("compat check (should match pairwise):")
for mat in [1,2]:
    for th in [0.3, 1.0, 2.5, -0.7, -2.0]:
        print(mat, th, check_disp(mat,1.0,th))
