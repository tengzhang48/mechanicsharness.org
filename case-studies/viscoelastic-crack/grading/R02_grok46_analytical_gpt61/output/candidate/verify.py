#!/usr/bin/env python3
"""Check the saved evidence without requiring a reference solution.

Run after solve.py: python verify.py case.json results.json
These are consistency/convergence checks, not claims of external accuracy.
"""
import json
import math
import sys
from pathlib import Path


def finite_tree(value):
    if isinstance(value, dict):
        return all(finite_tree(v) for v in value.values())
    if isinstance(value, list):
        return all(finite_tree(v) for v in value)
    if isinstance(value, float):
        return math.isfinite(value)
    return True


def main():
    inp = json.loads(Path(sys.argv[1] if len(sys.argv)>1 else 'case.json').read_text())
    result = json.loads(Path(sys.argv[2] if len(sys.argv)>2 else 'results.json').read_text())
    assert finite_tree(result), 'Non-finite JSON value'
    assert result['unresolved'] == []
    assert len(result['cases']) == len(inp['speeds'])
    assert result['endpoint']['speed_side_checks'][0]['limiting_opening_minus_target'] > 0
    assert result['endpoint']['speed_side_checks'][1]['limiting_opening_minus_target'] < 0
    assert result['endpoint']['independent_series_checks'][-1]['absolute_residual'] < 1e-9
    assert all(q['opening_shortfall'] > 0 for q in
               result['endpoint']['finite_length_sequence_at_Vc'])
    expected = {(float(x),float(y)) for x in inp['sample_X']
                for y in inp['sample_y_over_h']}
    for case, speed in zip(result['cases'],inp['speeds']):
        assert case['V'] == speed
        assert case['admissible'] and speed < result['critical_speed']
        assert 0 < case['a'] < 1 and case['epsilon'] > 0 and case['R'] > 1
        assert len(case['fields']) == len(expected)
        assert {(r['X'],r['y_over_h']) for r in case['fields']} == expected
        for r in case['fields']:
            assert all(r[k] is not None for k in ['T1','T2','E1','E2'])
            if r['y_over_h'] == 1:
                assert abs(r['T1']) < 1e-12 and abs(r['E1']) < 1e-12
        D = case['bulk_dissipation_ratio']
        balance = abs(D+1-case['R'])/max(1,case['R'])
        assert D > 0 and balance < 5e-5
        integration = case['bulk_integration']
        assert integration['last_refinement_absolute_change'] < 2e-5
        assert integration['maximum_tail_check_absolute_change'] < 1e-6
        assert integration['refinements'][-1]['filter_maximum_component_error'] < 1e-8
        assert all(q['raw_integral_core'] > 0 and q['raw_integral_tail'] >= 0
                   for q in integration['refinements'])
        verify = case['verification']
        for key in ['equilibrium_finite_difference',
                    'strain_compatibility_finite_difference',
                    'material_law_finite_difference']:
            seq = verify[key]
            assert seq[-1]['max_absolute_residual'] < 1e-6
            for p,q in zip(seq[:-1],seq[1:]):
                assert q['max_absolute_residual'] < .3*p['max_absolute_residual']
        assert verify['critical_delta_absolute_error'] < 1e-8
        assert verify['intact_U_max_abs'] < 1e-9
        assert verify['cohesive_traction_max_error'] < 1e-12
        assert verify['separated_traction_max_error'] < 1e-12
        assert verify['top_E1_max_abs'] < 1e-12
        assert abs(verify['top_H_minus_pi_a']) < 1e-12
        assert verify['far_ahead_T_error'] < 1e-9
        assert verify['far_ahead_E_error'] < 1e-9
        assert verify['far_behind_T_max_abs'] < 1e-8
        assert verify['far_behind_E_max_abs'] < 1e-8
        assert verify['sample_field_max_alternative_kernel_error'] < 1e-8
        assert verify['local_power_minus_free_energy_rate_vs_square_loss'][-1]['max_absolute_residual'] < 1e-6
        print(f'PASS V={speed:g}: 15 field rows; relative balance={balance:.3e}; '
              f'last bulk refinement={integration["last_refinement_absolute_change"]:.3e}; '
              f'tail change={integration["maximum_tail_check_absolute_change"]:.3e}')
    print('PASS: finite JSON, all requested rows, endpoint side/series checks, '
          'PDE/material/local-energy residuals, boundaries and critical separation.')


if __name__ == '__main__':
    main()
