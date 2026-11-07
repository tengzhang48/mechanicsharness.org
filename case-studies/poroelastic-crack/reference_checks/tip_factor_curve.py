"""Drained-tip factor K1(0+)/sqrt(pi L) = n11/(1+eta_t) against the coupling eta_t.

Runs the archived Codex Wilson factorization (evidence/corpus-runs/.../codex/scratch/wilson_factor.py,
unchanged) at three grid counts and extrapolates quadratically, exactly as that script does.
Writes tip_factor_curve.json next to this file.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

here = Path(__file__).resolve().parent
source = here.parent / 'grading/review_B_K10-codex/output/scratch/wilson_factor.py'
spec = importlib.util.spec_from_file_location('wilson_factor', source)
wf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wf)

rows = []
for eta in [0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 1 / 3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    n11 = [wf.wilson(eta, count, 0.01)['n11'] for count in (16384, 65536, 262144)]
    n21 = [wf.wilson(eta, count, 0.01)['n21'] for count in (16384, 65536, 262144)]
    w = np.array([1 / 3, -2, 8 / 3])
    e11, e21 = float(w @ n11), float(w @ n21)
    rows.append({'eta_t': eta, 'n11': e11, 'n21': e21, 'K_tip': e11 / (1 + eta),
                 'energy_identity_residual': e11**2 + e21**2 - (1 + eta),
                 'heuristic': 1 / (1 + eta), 'elastic_energy_only': (1 + eta) ** -0.5})
    print(json.dumps(rows[-1]), flush=True)
(here / 'tip_factor_curve.json').write_text(json.dumps({'source': str(source.relative_to(here.parent)),
                                                       'scale': 0.01, 'counts': [16384, 65536, 262144],
                                                       'rows': rows}, indent=1) + '\n')
