import numpy as np

def mode1_cart_stress(K, x, y):
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    pref = K/np.sqrt(2*np.pi*r)
    c, s = np.cos(th/2), np.sin(th/2)
    c3, s3 = np.cos(3*th/2), np.sin(3*th/2)
    sxx = pref*c*(1 - s*s3)
    syy = pref*c*(1 + s*s3)
    sxy = pref*s*c*c3
    return sxx, syy, sxy

def J_circle(K, R, E, nu, plane='stress', k=200001):
    a = np.linspace(-np.pi, np.pi, k)
    x, y = R*np.cos(a), R*np.sin(a)
    th = a
    if plane == 'stress':
        kappa = (3-nu)/(1+nu)
        # strains plane stress
        nu_eff = nu
    else:
        kappa = 3-4*nu
    G = E/(2*(1+nu))
    sxx, syy, sxy = mode1_cart_stress(K, x, y)
    if plane == 'stress':
        exx = (sxx - nu*syy)/E
        eyy = (syy - nu*sxx)/E
    else:
        # plane strain: ezz=0, szz = nu*(sxx+syy)
        szz = nu*(sxx+syy)
        exx = (1+nu)/E*((1-nu)*sxx - nu*syy)
        eyy = (1+nu)/E*((1-nu)*syy - nu*sxx)
    exy = (1+nu)/E*sxy
    W = 0.5*(sxx*exx + syy*eyy + 2*sxy*exy)  # szz*ezz=0 for plane stress (szz=0); plane strain ezz=0
    # displacement field
    pref = K/(2*G)*np.sqrt(R/(2*np.pi))
    c, s = np.cos(th/2), np.sin(th/2)
    fac = kappa - np.cos(th)
    ur = pref*(c*fac*np.cos(th) + s*fac*np.sin(th))  # ux cos + uy sin  -> radial
    ut = pref*(-c*fac*np.sin(th) + s*fac*np.cos(th))
    ux = pref*c*fac
    uy = pref*s*fac
    Nx, Ny = np.cos(a), np.sin(a)
    tx = sxx*Nx + sxy*Ny
    ty = sxy*Nx + syy*Ny
    # dux/dx via analytic? use finite diff in x with small h
    h = R*1e-6
    sxx2, syy2, sxy2 = mode1_cart_stress(K, x+h, y)
    # recompute ux at x+h
    r2 = np.hypot(x+h, y); th2 = np.arctan2(y, x+h)
    pref2 = K/(2*G)*np.sqrt(r2/(2*np.pi))
    c2, s2 = np.cos(th2/2), np.sin(th2/2)
    fac2 = kappa - np.cos(th2)
    ux2 = pref2*c2*fac2
    dux = (ux2-ux)/h
    integrand = (W*Nx - (tx*dux + ty*0))  # need duy/dx too
    # duy/dx
    uy2 = pref2*s2*fac2
    duy = (uy2-uy)/h
    integrand = (W*Nx - (tx*dux + ty*duy))
    return np.trapezoid(integrand*R, a)

for nu in [0.0, 0.2, 0.3, 0.45]:
    E = 70000.0; K = 10.0
    Jps = J_circle(K, 1.0, E, nu, 'stress')
    Jpe = J_circle(K, 1.0, E, nu, 'strain')
    print(nu, 'plane stress J*E/K^2=', Jps*E/K**2, ' plane strain J*E/((1-nu^2)K^2)=', Jpe*E/((1-nu**2)*K**2))
