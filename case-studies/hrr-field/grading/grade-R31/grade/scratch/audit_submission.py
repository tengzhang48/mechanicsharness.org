"""Local grading diagnostics; leaves submitted inputs unchanged."""
from pathlib import Path
import importlib.util
import json
import runpy

import numpy as np
from scipy.integrate import simpson

base = Path(__file__).parent
submitted = runpy.run_path(str(base / 'submission/ANSWER.py'))
n, p, k = (submitted[name] for name in ('n', 'p', 'k'))
rhs, solution, scale = (submitted[name] for name in ('rhs', 'sol', 'sc'))
theta = np.linspace(0, np.pi, 20001)
raw = solution(theta)
phi, v, w, z = raw * scale
f4 = np.asarray(rhs(theta, raw)[3]) * scale
rr, tt, rt = k*phi+w, k*(k-1)*phi, -(k-1)*v
d = rr-tt
dp = k*(2-k)*v+z
dpp = k*(2-k)*w+f4
rtp, rtpp = -(k-1)*w, -(k-1)*z
q2 = 0.75*d*d+3*rt*rt
q2p = 1.5*d*dp+6*rt*rtp
q2pp = 1.5*(dp*dp+d*dpp)+6*(rtp*rtp+rt*rtpp)
m = (n-1)/2
q = q2**m
qp = q*m*q2p/q2
qpp = q*(m*q2pp/q2+m*(m-1)*(q2p/q2)**2)
er = 0.75*q*d
g = er/p
gp = 0.75*(qp*d+q*dp)/p
gpp = 0.75*(qpp*d+2*qp*dp+q*dpp)/p
ell = (gp-3*q*rt)/(p+1)  # submitted reconstruction
ellp = (gpp-3*(qp*rt+q*rtp))/(p+1)
uxr = p*g*np.cos(theta)-(gp-ell)*np.sin(theta)
uxt = p*ell*np.cos(theta)-(ellp+g)*np.sin(theta)
energy = n/(n+1)*q2**((n+1)/2)
work_rr = 2*simpson(rr*uxr, x=theta)
work_rt = 2*simpson(rt*uxt, x=theta)
wcos = 2*simpson(energy*np.cos(theta), x=theta)
j_drop = wcos-work_rr
j_with_rt = j_drop-work_rt

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

reference = load(base / 'reference/reference.py', 'grading_reference')
redo = load(base / 'submission/redo.py', 'grading_redo')
probe = np.array([0.7, -0.2, -0.4, 0.1])
redo_rhs_error = float(np.max(np.abs(
    np.asarray(redo.rhs(0.5, probe))-reference.ode(0.5, probe, n, k)
)))

result = {
    'number_errors': {
        'p': abs(0.1333333-0.133333333333),
        'H': abs(2.15092-1.893569725213),
        'C': abs(0.29178-0.457305783974),
    },
    'normalized_face_stresses': {'rr': float(rr[-1]), 'tt': float(tt[-1]), 'rt': float(rt[-1])},
    'normalized_phi_third_at_zero': float(z[0]),
    'normalized_u_theta_at_zero': float(ell[0]),
    'max_correct_physical_shear_residual': float(np.max(np.abs(gp+(p-1)*ell-3*q*rt))),
    'max_hoop_strain_residual': float(np.max(np.abs(g+ellp+er))),
    'contour_decomposition_using_upper_half_and_mode_I_parity': {
        'energy_cos': float(wcos), 'rr_work': float(work_rr), 'rt_work': float(work_rt),
        'J_omitting_rt': float(j_drop), 'J_including_rt': float(j_with_rt),
        'note': 'Neither integral makes the submitted incompatible, asymmetric angular solution valid.',
    },
    'corrected_angular_setup_in_redo_max_rhs_difference_at_probe': redo_rhs_error,
    'energy_identity_check': {
        'max_W_minus_n_over_n_plus_1_Q_power': float(np.max(np.abs(energy-n/(n+1)*q2**((n+1)/2)))),
        'max_W_minus_Q_power': float(np.max(np.abs(energy-q2**((n+1)/2)))),
    },
}
(base / 'audit_results.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
