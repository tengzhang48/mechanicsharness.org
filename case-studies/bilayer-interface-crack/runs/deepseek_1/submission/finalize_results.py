"""
Final results for the bilayer interface-crack problem.

eps, G and |K| are evaluated in closed form (beam energy + the standard
bimaterial G-K relation).  The phase angles are obtained from the
high-order FEM crack-opening-displacement fits stored in data/cod_*.npz;
see method.md for how those files were generated and for the independent
checks (stress extrapolation, interaction integral, mesh convergence).
"""
import json, math
import numpy as np

# ---------------- material and geometry ----------------
nu = 1/3
E1, E2 = 9.0, 1.0
h, H = 1.0, 2.0
mu1 = E1/(2*(1+nu)); mu2 = E2/(2*(1+nu))
kappa = 3-4*nu
E1s = E1/(1-nu**2); E2s = E2/(1-nu**2)
Estr = 1.0/(0.5*(1/E1s + 1/E2s))

# ---------------- bimaterial constants ----------------
eps = math.log((kappa/mu1 + 1/mu2)/(kappa/mu2 + 1/mu1))/(2*math.pi)
beta = (mu1*(kappa-1)-mu2*(kappa-1))/(mu1*(kappa+1)+mu2*(kappa+1))

# ---------------- energy release rates (beam/J integral) ----------------
# G_A / (P^2/(E2 h)) and G_B / (M^2/(E2 h^3))
GA = (1.0/(2*E1s*h) + 1.0/(2*E2s*H)
      + 6.0*((h+H)/2.0)**2/(E2s*H**3))
GB = 6.0/(E1s*h**3) + 6.0/(E2s*H**3)
# with h=1, H=2 these are 287/162 and 34/27

# ---------------- |K| from G = (1-beta^2)|K|^2/E* ----------------
KA = math.sqrt(GA*Estr/(1-beta**2))
KB = math.sqrt(GB*Estr/(1-beta**2))

# ---------------- phase angles from FEM COD fits ----------------
# effective COD constant, including the 1/(1+2 i eps) factor
Ceff = 8.0/(Estr*math.cosh(math.pi*eps))/(1+2j*eps)

def fit_K(tag, rmin, rmax, deg):
    z = np.load(f'data/cod_{tag}.npz')
    r, d = z['r'], z['d']
    Q = d*np.sqrt(2*np.pi/r)/Ceff*r**(-1j*eps)
    sel = (r >= rmin) & (r <= rmax)
    rr, qq = r[sel], Q[sel]
    A = np.vander(rr, deg+1, increasing=True)
    cr = np.linalg.lstsq(A, qq.real, rcond=None)[0][0]
    ci = np.linalg.lstsq(A, qq.imag, rcond=None)[0][0]
    return cr+1j*ci

def phase(tag):
    # average of the two central polynomial fits (r in [0.01,0.1] and [0.02,0.1])
    K1 = fit_K(tag, 0.01, 0.1, 2)
    K2 = fit_K(tag, 0.02, 0.1, 2)
    K = 0.5*(K1+K2)
    return K, math.degrees(np.angle(K)), abs(K)

KA_fit, psiA, absA_fit = phase('A')
KB_fit, psiB, absB_fit = phase('B')
KC_fit, psiC, absC_fit = phase('C')

# The phase values we report (rounded to 0.01 deg; the spread between
# fit windows, meshes and independent extraction methods is < 0.15 deg).
psiA_r = round(psiA, 2)
psiB_r = round(psiB, 2)
psiC_r = round(psiC, 2)

# superposition check for M = P h
Ksup = KA*math.cos(math.radians(psiA_r)) + KB*math.cos(math.radians(psiB_r)) \
     + 1j*(KA*math.sin(math.radians(psiA_r)) + KB*math.sin(math.radians(psiB_r)))
psi_sup = math.degrees(np.angle(Ksup))

results = {
    "eps": eps,
    "caseA": {
        "G_E2_h_over_P2": GA,
        "Kabs_sqrt_h_over_P": KA,
        "psi_deg": psiA_r,
    },
    "caseB": {
        "G_E2_h3_over_M2": GB,
        "Kabs_h32_over_M": KB,
        "psi_deg": psiB_r,
    },
    "psi_both_deg": psiC_r,
}

with open('results.json', 'w') as f:
    json.dump(results, f, indent=2)
    f.write('\n')

if __name__ == '__main__':
    print('eps       =', eps)
    print('beta      =', beta, ' E* =', Estr)
    print('G_A       =', GA, '   G_B =', GB)
    print('|K_A|     =', KA, '   |K_B| =', KB)
    print('COD fits:')
    print('  A: psi = %.4f deg, |K|_fit = %.5f (exact %.5f)' % (psiA, absA_fit, KA))
    print('  B: psi = %.4f deg, |K|_fit = %.5f (exact %.5f)' % (psiB, absB_fit, KB))
    print('  C: psi = %.4f deg, |K|_fit = %.5f' % (psiC, absC_fit))
    print('superposition A+B -> psi = %.4f deg (direct combined fit %.4f)' % (psi_sup, psiC))
    print('reported: psi_A=%.2f, psi_B=%.2f, psi_both=%.2f' % (psiA_r, psiB_r, psiC_r))
    print(open('results.json').read())
