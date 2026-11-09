"""Read-only comparison of the copied submission, diagnostic repair, and key."""
import contextlib
import io
import json
import runpy
from pathlib import Path
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    submitted = runpy.run_path('scratch/solve.py')
    corrected = runpy.run_path('scratch/solve_even_derivative.py')

reference = json.loads(Path('reference_output.json').read_text())['route_A_angular_ODE']
answer = {
    'slip_0.02_mm': 1.83e-3,
    'slip_0.08_mm': 2.49e-3,
    'tau_rz_0.02mm_60deg_MPa': 42.4,
    'tau_thetaz_0.02mm_60deg_MPa': 151.7,
    'gamma_contour_0deg_mm': .07167,
    'gamma_contour_90deg_mm': .03831,
    'gamma_contour_180deg_mm': .02048,
}

def evaluated_values(namespace):
    vals = {
        'slip_0.02_mm': namespace['Dw1'],
        'slip_0.08_mm': namespace['Dw2'],
        'tau_rz_0.02mm_60deg_MPa': namespace['trz'],
        'tau_thetaz_0.02mm_60deg_MPa': namespace['ttz'],
    }
    for deg in (0, 90, 180):
        q = float(namespace['Qval'](np.deg2rad(deg)))
        vals[f'gamma_contour_{deg}deg_mm'] = (namespace['A']*q/.01)**(1/(1-namespace['beta']))
    return vals

submitted_values = evaluated_values(submitted)
corrected_values = evaluated_values(corrected)
theta = -np.pi/3
h = 1e-6
numerical_derivative = float((submitted['fval'](theta+h)-submitted['fval'](theta-h))/(2*h))
checks = {
    'requested_values': {
        key: {
            'answer': answer[key],
            'code': float(submitted_values[key]),
            'reference': reference[key],
            'answer_relative_error_percent': 100*abs(answer[key]/reference[key]-1),
            'code_relative_error_percent': 100*abs(submitted_values[key]/reference[key]-1),
            'within_two_percent': bool(abs(answer[key]/reference[key]-1) <= .02),
            'corrected_parity_relative_error': float(abs(corrected_values[key]/reference[key]-1)),
        } for key in reference
    },
    'negative_sixty_degree_derivative': {
        'finite_difference_of_fval': numerical_derivative,
        'submitted_uval': float(submitted['uval'](theta)),
        'correct_even_uval': float(corrected['uval'](theta)),
    },
    'flux_coefficient_submitted': float(submitted['alpha']),
    'flux_coefficient_corrected': float(corrected['alpha']),
    'amplitude_submitted': float(submitted['A']),
    'amplitude_corrected': float(corrected['A']),
    'physical_flux_using_submitted_amplitude': float(corrected['J_of'](submitted['A'], .02)),
    'face_derivative_residual': float(submitted['upi']),
    'max_corrected_parity_relative_error': float(max(abs(corrected_values[k]/reference[k]-1) for k in reference)),
}
Path('scratch/grading_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
