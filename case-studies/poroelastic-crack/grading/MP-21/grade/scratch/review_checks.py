"""Independent grading checks; do not modify or import the submission."""
import json
from pathlib import Path

import numpy as np

nu, nu_u, alpha = 0.2, 0.4, 0.89
G = L = P = 1.0
beta = alpha * (1 - 2 * nu) / (2 * G)
eta_back = alpha * (1 - 2 * nu) / (2 * (1 - nu))
key = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": np.pi * (1 - nu_u),
    "V1_inf": np.pi * (1 - nu),
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (nu_u - nu) / (alpha * (1 - 2 * nu)),
    "K2_0": 0.0,
    "K2_inf": -eta_back,
    "V2_inf": -np.pi * alpha * (1 - 2 * nu) / 2,
}
reported = json.loads(Path('submission/results.json').read_text())
print('Reported values against corrected key (1% nonzero, 0.01 for zero):')
for name, expected in key.items():
    value = reported[name]
    tol = 0.01 * abs(expected) if expected else 0.01
    print(f'{name:16s} submitted={value:.10f} key={expected:.10f} '
          f'pass={abs(value - expected) <= tol}')

print('\nExact swelling-disk counterexample to p -> 1 implies sigma -> 0:')
print('p=1 for r<R, p=0 for r>R; u=a*r inside, u=a*R^2/r outside.')
# Obtained by displacement and radial-traction continuity at r=R.
a = alpha * (1 - 2 * nu) / (4 * G * (1 - nu))
stress_inside = 2 * G * a / (1 - 2 * nu) - alpha
for radius in (2.0, 20.0, 200.0):
    boundary_stress_outside = -2 * G * a
    displacement_jump = a * radius - a * radius**2 / radius
    print(f'R={radius:6.1f}  p(r=1)=1  sigma_rr(r=1)={stress_inside:.8f} '
          f'traction_jump={stress_inside - boundary_stress_outside:.3e} '
          f'displacement_jump={displacement_jump:.3e}')
print(f'strain inside={a:.8f}; free swelling strain claimed={beta:.8f}')
print('Every finite disk solution has decaying stresses and displacement at infinity.')
print('Taking R to infinity pointwise leaves compressive stress -eta, not zero.')
print(f'A traction-free crack in this prestress requires face stress +eta: '
      f'K2_inf={-eta_back:.8f}, V2_inf={-np.pi*(1-nu)*eta_back:.8f}')

print('\nCrack-volume reciprocity check for expanding pressure plateaus:')
print('Use the unit-pressure drained Griffith stress trace S_aux.')
print('For zero mechanical face traction, V = beta*integral(p*S_aux dA).')
print('In elliptic coordinates z=L*cosh(mu+i*theta), '
      'integral_0^(2*pi)(S_aux*Jacobian)dtheta = -2*pi*L^2*exp(-2*mu).')
print('For p=1 inside mu<M and p=0 outside, the exact finite-domain integral is:')
for M in (1.0, 2.0, 4.0, 8.0):
    volume = -beta*np.pi*L**2*(1-np.exp(-2*M))
    print(f'M={M:.1f}, V={volume:.10f}')
print(f'The uniform pointwise pressure limit leaves V={-beta*np.pi*L**2:.10f}, not zero.')

print('\nCoupled pressure equation from compatibility and mass balance:')
B = alpha * (1 - 2 * nu) / (1 - nu)
forcing = (nu_u - nu) * (1 - nu) / (alpha * (1 - 2 * nu) * (1 - nu_u))
print(f'h = sigma_kk + {B:.8f}*p is harmonic; '
      f'p_t = c*laplacian(p) - {forcing:.8f}*h_t.')
print('The submission check 3b sets h_t=0 without solving the mechanical crack problem.')

print('\nWestergaard branch audit using the submitted formula:')
def submitted_stress(x, y):
    z = complex(x, y)
    s = np.sqrt(z*z - L*L + 0j)
    Z = P * (z/s - 1)
    dZ = -P * L**2 / (z*z - L*L + 0j)**1.5
    return np.array([Z.real-y*dZ.imag, Z.real+y*dZ.imag, -y*dZ.real])
for x in (2.0, -2.0, 100.0, -100.0):
    print(f'x={x:6.1f}, y=0.5: sigma_xx, sigma_yy, tau_xy={submitted_stress(x, 0.5)}')

print('\nSource of the DD-volume quadrature discrepancy:')
x = np.linspace(-L, L, 801)
ellipse = 2*(1-nu)*np.sqrt(np.maximum(L*L-x*x, 0))
integrate = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
print(f'801-point trapezoid of the exact opening={integrate(ellipse, x):.8f}; '
      f'exact volume={np.pi*(1-nu):.8f}')

print('\nShort-time coefficients:')
faces = 16*(1-nu_u)*(nu_u-nu)/(np.sqrt(np.pi)*(1-nu))
print(f'submitted local coefficient={reported["V1_short_coeff"]:.9f}; '
      f'key face contribution={faces:.9f}; corrected full coefficient={key["V1_short_coeff"]:.9f}')
eta_u = (nu_u-nu)/(alpha*(1-2*nu))
kappaG_over_c = alpha**2*(1-2*nu)**2*(1-nu_u)/(2*(1-nu)*(nu_u-nu))
reciprocal_faces = 32*eta_u**2*kappaG_over_c/np.sqrt(np.pi)
print(f'Independent face-only reciprocal calculation: '
      f'A_faces = 32*eta_u^2*kappa*G/(c*sqrt(pi)) = {reciprocal_faces:.9f}.')

print('\nHydraulic loading 1+2 endpoint values using the key:')
print(f'K(0+)={key["K1_0"]+key["K2_0"]:.8f}; '
      f'K(inf)={key["K1_inf"]+key["K2_inf"]:.8f}; '
      f'V(inf)={key["V1_inf"]+key["V2_inf"]:.8f}')
