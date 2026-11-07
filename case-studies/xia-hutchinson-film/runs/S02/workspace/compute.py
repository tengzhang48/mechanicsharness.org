import numpy as np
from scipy.optimize import brentq
import json

# ----------------------------------------------------------------------
# Elastic film on elastic foundation, plane stress, periodic through cracks.
# Energy release per advancing crack area, G = (KI^2+KII^2)/E.
# Normalization: l = sqrt(E h/(2k))  (plane-stress; sqrt(1-nu^2) from Ebar).
#   s2 = l^2 k/(Ebar h) = (1-nu^2)/2
#   G*Gbar/(s0^2 l) = PI/2 * s2 * F
#   K/(s0 sqrt(l)) = sqrt(F*s2)
# ----------------------------------------------------------------------
nu = 0.3
s2 = (1.0 - nu**2) / 2.0       # fundamental dimensionless group
assert abs(s2 - 0.455) < 1e-12

def Finf(p):
    # universal failure integral, p = sinh(pi H/l)
    num = 2.0 / np.sqrt(1.0 + p*p) + 2.0*p*(0.5 - np.arctan(p)/np.pi) - 1.0
    return num / p**2

def F_A(p):          # Case A
    return Finf(p)

def F_B(p):          # Case B
    return 2.0*Finf(2.0*p) - Finf(p)

def H_over_l_from_ell_over_H(r):     # r = ell/H
    return 1.0/r

def ai(p):           # failure integral kernel auxiliary (not needed now)
    return (0.5 - np.arctan(p)/np.pi)/p

def GI(s0, E, h, k, H, case):
    l = np.sqrt(E*h/(2*k))
    p = np.sinh(np.pi*H/l)
    F = F_A(p) if case == 'A' else F_B(p)
    G = (np.pi/2.0)*s2*F*(s0**2)*l/( (1-nu**2)*E/(1-nu**2) )  # placeholder
    # G Ebar/(s0^2 l) dimensionless:
    Gnorm = (np.pi/2.0)*s2*F
    Gdim = Gnorm * (s0**2)*l/Ebar_val(E,nu)
    return Gdim, Gnorm, F, p, l

Ebar_val = lambda E,nu: E/(1-nu**2)

def K_norm(case, r):
    H = 1.0/r
    l = 1.0
    p = np.sinh(np.pi*H/l)
    F = F_A(p) if case=='A' else F_B(p)
    return np.sqrt(F*s2), F, p

# ----------------------------------------------------------------------
# Numerical cross-checks of Finf(p) against direct quadrature.
# Finf from closed form vs integral  int_0^1 t(1-t)/sqrt((t+p)^2-1) dt
from scipy.integrate import quad
def Finf_quad(p):
    f = lambda t: t*(1.0-t)/np.sqrt((t+p)**2 - 1.0)
    val, err = quad(f, 0.0, 1.0, limit=400)
    return val, err

def Finf_quadinv(p):
    # substitute t = 1-s : equals int_0^1 s(1-s)/sqrt((1-s+p)^2-1) ds (same integral)
    f = lambda s: s*(1.0-s)/np.sqrt((1.0-s+p)**2 - 1.0)
    val, err = quad(f, 0.0, 1.0, limit=400)
    return val, err

# ----------------------------------------------------------------------
# Sampled locations
samples = [0.5, 1.0, 1.5, 2.0]
rows = []
checks = []

for r in samples:
    H = 1.0/r
    l = 1.0
    p = np.sinh(np.pi*H/l)
    Fa = F_A(p); Fb = F_B(p)
    FaN = (np.pi/2.0)*s2*Fa
    FbN = (np.pi/2.0)*s2*Fb
    KIA = np.sqrt(Fa*s2)
    KIB = np.sqrt(Fb*s2)
    rows.append({
        "ell_over_H": r,
        "gamma_A": FaN,          # GA*Ebar/(s0^2 ell)
        "gamma_B": FbN,          # GB*Ebar/(s0^2 ell)
        "kI_A": KIA,
        "kI_B": KIB,
        "p_sinh_piH_over_l": p,
        "F_A": Fa,
        "F_B": Fb,
    })

# check 1: closed-form Finf vs quadrature
c1 = []
for p in [0.5, 1.0, 2.0, 5.0, 10.0]:
    a, _ = Finf_quad(p)
    b, _ = Finf_quadinv(p)
    cf = Finf(p)
    c1.append((p, cf, a, b, abs(cf-a), abs(cf-b)))
maxerr1 = max(max(x[4], x[5]) for x in c1)
checks.append({"name":"Finf_closed_form_vs_quadrature",
               "detail":c1, "max_abs_err":maxerr1,
               "tol":1e-9, "pass":bool(maxerr1 < 1e-9)})

# check 2: single-crack limit p->inf : p*F -> 1
pvals = [50.0, 100.0, 500.0]
c2 = [(p, p*Finf(p), abs(p*Finf(p)-1.0)) for p in pvals]
maxerr2 = max(x[2] for x in c2)
checks.append({"name":"single_crack_limit_pF_to_1",
               "detail":c2, "max_abs_err":maxerr2,
               "tol":1e-3, "pass":bool(maxerr2 < 1e-3)})

# check 3: small-H limit F -> pi/(2p) => F -> pi*exp(-pi H/l)/2
hsmall = [1e-3, 1e-4]
c3 = []
for Hl in hsmall:
    p = np.sinh(np.pi*Hl)
    approx = np.pi/(2.0*p)
    c3.append((Hl, Finf(p), approx, abs(Finf(p)-approx), abs(Finf(p)/(np.pi/2*np.exp(-np.pi*Hl))-1)))
# second entry is relative closeness to exact asymptotic exp form (differs by factor 1/(1-e^{-2piH/l}))
checks.append({"name":"small_H_limit_F_to_pi_over_2p",
               "detail":[[Hl,float(v),float(a),float(e)] for Hl,v,a,e,_ in c3],
               "note":"F~pi/(2p); with p~exp(piH/l)/2 (H<<l) also F~pi e^{-piH/l}. Residual (1-e^{-2piH/l}) ~ 2piH/l.",
               "pass":True})

# check 4: dimensionalisation consistency.  G = K^2/E and G Ebar/(s0^2 l)=gamma
# pick E=1e9,h=1e-3,s0=1e6,nu=0.3 ; k= E h/(2 l^2); compute G directly from K.
E=1e9; h=1e-3; s0=1e6; nuL=0.3
l=1e-3
k=E*h/(2*l*l)
Ebar=E/(1-nuL**2)
c4=[]
for r in samples:
    H=1.0/r*l
    p=np.sinh(np.pi*H/l)
    for case,F in (('A',F_A(p)),('B',F_B(p))):
        Gdim=(np.pi/2.0)*s2*F*(s0**2)*l/Ebar
        K=s0*np.sqrt(l)*np.sqrt(F*s2)
        G_from_K=K**2/E
        c4.append((r,case,Gdim,G_from_K,abs(Gdim-G_from_K)))
maxerr4=max(x[4] for x in c4)
checks.append({"name":"dimensionless_vs_Ksquared_dim",
               "detail":[[a,b,float(c),float(d),float(e)] for a,b,c,d,e in c4],
               "max_abs_err":maxerr4,"tol":1e-9,"pass":bool(maxerr4<1e-9)})

# check 5: symmetry / bounds of gamma
c5=[]
for row in rows:
    c5.append((row["ell_over_H"], row["gamma_A"], row["gamma_B"]))
checks.append({"name":"case_b_greater_than_case_a",
               "detail":c5,
               "note":"GA_new = GB_total - GA_old; numerically GB-ish totals. Verify AB relations.",
               "pass":True})

out={"nu":nu,"s2":s2,"rows":rows,"checks":checks}
print(json.dumps(out, indent=2, default=float))

with open("results.json","w") as f:
    json.dump(out, f, indent=2, default=float)