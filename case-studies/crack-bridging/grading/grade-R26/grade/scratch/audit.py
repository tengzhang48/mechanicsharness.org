"""Grader audit of the submitted event and supplementary diagnostics."""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import solve_lib as sl

key = json.loads(Path('../reference_results.json').read_text())['answer']
Le = brentq(lambda L: sl.delta_end(L, sl.Mstar) - sl.dc,
            0.010, 0.011, xtol=1e-16, rtol=8.9e-16)
xx = np.linspace(0.0, Le, 5001)
delta, coefficients = sl.solution(Le, sl.Mstar, xx)
P, dP, d2P, d3P = sl.basis(xx)
curvature = coefficients @ d2P
slope = coefficients @ dP
third = coefficients @ d3P
tip_expected = np.sqrt(2 * sl.b * sl.CM * sl.G0)
invariant = (curvature**2 / (2 * sl.b * sl.CM)
             - slope * third / (sl.b * sl.CM)
             - sl.G0 - sl.bridging_work(delta))
computed = {
    'M_each_arm_N_m': float(sl.Mstar),
    'L_mm': float(Le * 1e3),
    'delta_at_L_over_2_um': float(sl.solution(Le, sl.Mstar, [Le/2])[0][0] * 1e6),
}
reported = {
    'main_bold_results': dict(zip(key, [0.2126, 10.80, 38.2])),
    'parenthetical_results': dict(zip(key, [0.212626, 10.8001, 38.219])),
    'unrounded_code_results': computed,
}

def roots_at(L):
    d0 = sl.delta_end(L, 0.0)
    q = sl.delta_end(L, 1.0) - d0
    polynomial = [sl.CM/(2*sl.b) + sl.sig*q*q/(2*sl.dc),
                  -q*sl.sig*(1-d0/sl.dc),
                  -(sl.G0+sl.bridging_work(d0))]
    return sorted(float(r.real) for r in np.roots(polynomial)
                  if abs(r.imag) < 1e-10 and r.real > 0)

branches = []
for M in roots_at(Le):
    d, c = sl.solution(Le, M, xx)
    branches.append({'M_N_m': M, 'delta_rear_um': float(d[-1] * 1e6),
                     'tip_curvature_per_m': float((c @ d2P)[0]),
                     'minimum_gap_um': float(d.min() * 1e6)})

continuation = []
previous = sl.Minit
for L in np.linspace(1e-4, Le, 1200):
    M = min(roots_at(L), key=lambda v: abs(v-previous))
    rear = sl.delta_end(L, M)
    q = sl.delta_end(L, 1.0)-sl.delta_end(L, 0.0)
    derivative_M = sl.CM*M/sl.b-q*sl.sig*(1-rear/sl.dc)
    epsilon = 1e-8
    derivative_delta_L = (sl.delta_end(L+epsilon, M)
                          - sl.delta_end(L-epsilon, M))/(2*epsilon)
    derivative_L = -sl.sig*(1-rear/sl.dc)*derivative_delta_L
    continuation.append([L, M, rear, derivative_M, derivative_L])
    previous = M

off_L, off_M = 1.3*Le, 1.1*sl.Mstar
off_gap = sl.solution(off_L, off_M, np.linspace(0, off_L, 5001))[0]
off_traction = sl.sig*(1-off_gap/sl.dc)
result = {
    'reference': key,
    'reported_and_rerun': reported,
    'absolute_relative_errors_percent': {
        label: {name: abs(value/key[name]-1)*100 for name, value in values.items()}
        for label, values in reported.items()
    },
    'event_checks': {
        'actual_tip_curvature_per_m': float(curvature[0]),
        'expected_positive_tip_curvature_per_m': float(tip_expected),
        'rear_curvature_per_m': float(curvature[-1]),
        'expected_rear_curvature_per_m': float(sl.CM*sl.Mstar),
        'max_first_integral_residual_J_m2': float(np.max(np.abs(invariant))),
        'profile_nondecreasing': bool(np.all(np.diff(delta) > 0)),
        'delta_rear_um': float(delta[-1]*1e6),
        'delta_min_um': float(delta.min()*1e6),
    },
    'supplementary_M_of_L_roots_at_event_length': branches,
    'pre_event_continuation': {
        'M_increases': bool(np.all(np.diff(np.array(continuation)[:,1]) > 0)),
        'rear_gap_increases': bool(np.all(np.diff(np.array(continuation)[:,2]) > 0)),
        'minimum_partial_F_partial_M': float(min(row[3] for row in continuation)),
        'maximum_partial_F_partial_L_before_event': float(max(row[4] for row in continuation[:-1])),
        'endpoint': continuation[-1],
    },
    'supplementary_off_event_energy_sample': {
        'L_mm': float(off_L*1e3), 'M_N_m': float(off_M),
        'minimum_gap_um': float(off_gap.min()*1e6),
        'maximum_gap_um': float(off_gap.max()*1e6),
        'traction_min_MPa': float(off_traction.min()*1e-6),
        'traction_max_MPa': float(off_traction.max()*1e-6),
    },
    'bending_only_same_M_and_L_midpoint_um': float(sl.CM*sl.Mstar*(Le/2)**2/2*1e6),
    'bending_only_same_rear_opening_midpoint_um': sl.dc/4*1e6,
}
print(json.dumps(result, indent=2))
