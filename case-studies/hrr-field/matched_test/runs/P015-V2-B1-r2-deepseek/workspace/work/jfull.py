import numpy as np
E, nu, sy, J = 70000.0, 0.30, 10.0, 0.800
T = 4816.472187302733
beta = 1.0 / T - 1.0 / E
G = T / (2 * (1 + (0.5 - 0.2 * T / E)))
ne = 0.5 - 0.2 * T / E
A_hyd = (1 - 2 * nu) / (3 * E)
Ee = T
KA = np.sqrt(J * Ee / (2 * np.pi))
kappa = (3 - ne) / (1 + ne)

def polar(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    srr = KA*(1.25*np.cos(t/2)-0.25*np.cos(1.5*t))/np.sqrt(r)
    stt = KA*(0.75*np.cos(t/2)+0.25*np.cos(1.5*t))/np.sqrt(r)
    srt = KA*(0.25*np.sin(t/2)+0.25*np.sin(1.5*t))/np.sqrt(r)
    return r, t, srr, stt, srt

def stress_cart(x, y):
    r, t, srr, stt, srt = polar(x, y)
    c, s = np.cos(t), np.sin(t)
    sxx = srr*c*c + stt*s*s - 2*srt*s*c
    syy = srr*s*s + stt*c*c + 2*srt*s*c
    sxy = (srr-stt)*s*c + srt*(c*c-s*s)
    return sxx, syy, sxy

def strain_cart_full(x, y):
    r, t, srr, stt, srt = polar(x, y)
    S = srr + stt
    szz = -S/3
    m = (srr+stt+szz)/3
    d = [srr-m, stt-m, srt, szz-m]
    se = np.sqrt(1.5*(d[0]**2 + d[1]**2 + d[3]**2 + 2*d[2]**2))
    f = np.maximum(1 - sy/se, 0.0)
    err = (1+nu)/E*srr - nu/E*S + 1.5*beta*f*d[0]
    ett = (1+nu)/E*stt - nu/E*S + 1.5*beta*f*d[1]
    ert = (1+nu)/E*srt + 1.5*beta*f*d[2]
    exx = err*np.cos(t)**2 + ett*np.sin(t)**2 - 2*ert*np.sin(t)*np.cos(t)
    eyy = err*np.sin(t)**2 + ett*np.cos(t)**2 + 2*ert*np.sin(t)*np.cos(t)
    exy = (err-ett)*np.sin(t)*np.cos(t) + ert*(np.cos(t)**2-np.sin(t)**2)
    return exx, eyy, exy

def W_full(x, y):
    r, t, srr, stt, srt = polar(x, y)
    S = srr + stt
    szz = -S/3
    m = S/3
    d = [srr-m, stt-m, srt, szz-m]
    se = np.sqrt(1.5*(d[0]**2 + d[1]**2 + d[3]**2 + 2*d[2]**2))
    W = 0.5*((1+nu)/E*(srr**2+stt**2+szz**2+2*srt**2) - nu/E*S**2)
    W += 0.5*beta*np.maximum(se - sy, 0.0)**2
    return W

def u_cart(x, y):
    r = np.hypot(x, y); t = np.arctan2(y, x)
    ur = (KA/(2*G))*np.sqrt(r)*(kappa-np.cos(t))*np.cos(t/2)
    ut = -(KA/(2*G))*np.sqrt(r)*(kappa-np.cos(t))*np.sin(t/2)
    return ur*np.cos(t)-ut*np.sin(t), ur*np.sin(t)+ut*np.cos(t)

def J_circle(R, k=400001):
    a = np.linspace(-np.pi+1e-6, np.pi-1e-6, k)
    x, y = R*np.cos(a), R*np.sin(a)
    sxx, syy, sxy = stress_cart(x, y)
    W = W_full(x, y)
    h = 1e-5*R
    u1 = np.array(u_cart(x+h, y)); u0 = np.array(u_cart(x-h, y))
    dux = (u1[0]-u0[0])/(2*h); duy = (u1[1]-u0[1])/(2*h)
    Nx, Ny = np.cos(a), np.sin(a)
    tx, ty = sxx*Nx+sxy*Ny, sxy*Nx+syy*Ny
    return np.trapezoid((W*Nx-(tx*dux+ty*duy))*R, a)

print("Full-law J (leading stress + effective-leading displacements):")
for R in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2):
    print("  R=%.4f  J=%.8f" % (R, J_circle(R)))
