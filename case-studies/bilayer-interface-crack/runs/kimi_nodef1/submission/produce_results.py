#!/usr/bin/env python3
"""Generate results.json from the exact energy balance and the final FE fit."""
import json
import math

import numpy as np

from extract_k import extrapolate, read_k

SOURCE = 'q9_120_L20_n128.json'
q, r, Ks = read_k(SOURCE)
eps = q['parameters']['epsilon']

fitted = []
for K in Ks:
    K0, max_error, count = extrapolate(r, K, eps, rmax=0.02,
                                      model=2, rmin=0.001)
    fitted.append(K0)

# Exact energy-release rates from the far-behind beam energies.
GA = 287.0 / 162.0
GB = 34.0 / 27.0
# With beta=1/5 and D=640/(81 E2), G = (64/135) A^2/E2.
amplitude_A = math.sqrt(135.0 * GA / 64.0)
amplitude_B = math.sqrt(135.0 * GB / 64.0)

# The FE determines the direction of each complex SIF. Replace its remaining
# sub-percent amplitude discretization error by the exact energy amplitude.
phase_A = math.degrees(np.angle(fitted[0]))
phase_B = math.degrees(np.angle(fitted[1]))
KA = amplitude_A * np.exp(1j * math.radians(phase_A))
KB = amplitude_B * np.exp(1j * math.radians(phase_B))
phase_both = math.degrees(np.angle(KA + KB))
phase_decade = math.degrees(eps * math.log(10.0))

results = {
    'caseA': {
        'G_E2_h_over_P2': GA,
        'amplitude_sqrt_h_over_P': amplitude_A,
        'phase_at_h_deg': phase_A,
    },
    'caseB': {
        'G_E2_h3_over_M2': GB,
        'amplitude_h32_over_M': amplitude_B,
        'phase_at_h_deg': phase_B,
    },
    'phase_at_h_both_deg': phase_both,
    'phase_change_per_decade_deg': phase_decade,
}
with open('results.json', 'w') as f:
    json.dump(results, f, indent=2)
    f.write('\n')

print(json.dumps(results, indent=2))
print('FE/exact amplitude ratios before exact normalization:',
      abs(fitted[0]) / amplitude_A, abs(fitted[1]) / amplitude_B)
