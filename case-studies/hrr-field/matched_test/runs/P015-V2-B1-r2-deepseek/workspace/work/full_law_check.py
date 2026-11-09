import numpy as np
E, nu, sy, J = 70000.0, 0.30, 10.0, 0.800
T = 4816.472187302733
beta = 1.0 / T - 1.0 / E
ne = 0.5 - 0.2 * T / E
KA = np.sqrt(J * T / (2 * np.pi))

def full_strain(r, t):
    """Full nonlinear law applied to the leading stress field sigma = r^-1/2 Sigma(theta)."""
    c2, c32, s2, s32 = np.cos(t/2), np.cos(1.5*t), np.sin(t/2), np.sin(1.5*t)
    srr = KA*(1.25*c2-0.25*c32)/np.sqrt(r)
    stt = KA*(0.75*c2+0.25*c32)/np.sqrt(r)
    srt = KA*(0.25*s2+0.25*s32)/np.sqrt(r)
    S = srr+stt
    szz = 0.0
    m = S/3
    d0, d1, d2, d3 = srr-m, stt-m, srt, szz-m
    se = np.sqrt(1.5*(d0**2+d1**2+d3**2+2*d2**2))
    f = max(1-sy/se, 0.0)
    err = (1+nu)/E*srr - nu/E*S + 1.5*beta*f*d0
    ett = (1+nu)/E*stt - nu/E*S + 1.5*beta*f*d1
    ert = (1+nu)/E*srt + 1.5*beta*f*d2
    ezz = -nu/E*S + 1.5*beta*f*d3
    return err, ett, ert, ezz, se, f

print("theta = 0:  full-law sqrt(r)*eps_zz and bracket factor")
for r in (0.01, 0.001, 1e-4, 1e-5):
    err, ett, ert, ezz, se, f = full_strain(r, 0.0)
    print("  r=%.0e  sqrt(r)*ezz=%.9f  se=%.3f  bracket=%.6f" % (r, np.sqrt(r)*ezz, se, f))

t60 = np.pi/3
A_dev = 1.5/T + (nu-0.5)/E
A_hyd = (1-2*nu)/(3*E)

def eff_strain(t):
    c2, c32, s2, s32 = np.cos(t/2), np.cos(1.5*t), np.sin(t/2), np.sin(1.5*t)
    srr = (1.25*c2-0.25*c32); stt = (0.75*c2+0.25*c32); srt = (0.25*s2+0.25*s32)
    S = srr+stt
    err = A_hyd*S + A_dev*(srr-S/3)
    ett = A_hyd*S + A_dev*(stt-S/3)
    ert = A_dev*srt
    return KA*err, KA*ett, KA*2*ert

e_rr, e_tt, e_gr = eff_strain(t60)
print("\ntheta=60deg: effective-linear coefficients vs full-law limit")
print("  effective:  C_rr=%.9f C_tt=%.9f C_grth=%.9f" % (e_rr, e_tt, e_gr))
for r in (0.01, 0.001, 1e-4, 1e-5):
    err, ett, ert, ezz, se, f = full_strain(r, t60)
    # limit sqrt(r)*component; full law adds O(r^0)*sqrt(r) -> 0
    print("  r=%.0e  sqrt(r)*err=%.9f sqrt(r)*ett=%.9f sqrt(r)*2ert=%.9f bracket=%.6f"
          % (r, np.sqrt(r)*err, np.sqrt(r)*ett, np.sqrt(r)*2*ert, f))
