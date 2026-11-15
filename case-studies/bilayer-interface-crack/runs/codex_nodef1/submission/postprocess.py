"""Check the saved elasticity runs and write the requested results.json."""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from solve import BETA, EPSILON, G_COEFFICIENT, NU, exact_values


def read_run(name):
    return json.loads(Path(name).read_text())


def check_beam_energy(p, moment):
    """Independently integrate the stated remote stress energy densities."""
    q, w = leggauss(8)
    y1 = (q + 1) / 2
    y2 = q - 1
    s1 = -p - 12 * moment * (y1 - .5)
    s2 = p / 2 + 12 * (moment + 1.5 * p) * (y2 + 1) / 8
    return float((1 - NU**2) / 2 * (np.dot(w / 2, s1**2) / 9 + np.dot(w, s2**2)))


def main():
    ga, gb, gab = exact_values()
    amplitudes = np.sqrt(np.array([ga, gb]) / G_COEFFICIENT)
    fine = read_run('fem_fine.json')
    reference_key = 'average_0.0001_0.01'
    fitted = fine['fits'][reference_key]
    phase = np.asarray(fitted['phase'])
    # Use the exact energy-based amplitudes in superposing the two unit loads.
    # h=1, so M=P h corresponds to adding these two complex coefficients.
    kh = amplitudes * np.exp(1j * np.deg2rad(phase))
    both_phase = float(np.angle(kh.sum(), deg=True))

    results = {
        'caseA': {
            'G_E2_h_over_P2': ga,
            'amplitude_sqrt_h_over_P': float(amplitudes[0]),
            'phase_at_h_deg': round(float(phase[0]), 2),
        },
        'caseB': {
            'G_E2_h3_over_M2': gb,
            'amplitude_h32_over_M': float(amplitudes[1]),
            'phase_at_h_deg': round(float(phase[1]), 2),
        },
        'phase_at_h_both_deg': round(both_phase, 2),
        'phase_change_per_decade_deg': float(np.rad2deg(EPSILON * np.log(10))),
    }
    assert set(results) == {'caseA', 'caseB', 'phase_at_h_both_deg',
                            'phase_change_per_decade_deg'}
    required_fields = {
        'caseA': {'G_E2_h_over_P2', 'amplitude_sqrt_h_over_P', 'phase_at_h_deg'},
        'caseB': {'G_E2_h3_over_M2', 'amplitude_h32_over_M', 'phase_at_h_deg'},
    }
    for case, fields in required_fields.items():
        assert set(results[case]) == fields
        assert all(type(value) in (float, int) and np.isfinite(value)
                   for value in results[case].values())

    integrated_energies = [check_beam_energy(1, 0), check_beam_energy(0, 1),
                           check_beam_energy(1, 1)]
    analytic_energies = [ga, gb, ga + gb + 2 * gab]
    assert np.allclose(integrated_energies, analytic_energies, rtol=1e-13)

    # Polarization of G(P,M) determines the cosine of the phase difference.
    relative_phase_exact = float(np.rad2deg(np.arccos(gab / np.sqrt(ga * gb))))
    relative_phase_fem = float(phase[0] - phase[1])
    relative_phase_error = relative_phase_fem - relative_phase_exact
    assert abs(relative_phase_error) < .05

    amplitude_error = (np.asarray(fitted['amplitude']) / amplitudes - 1) * 100
    assert np.max(np.abs(amplitude_error)) < .1

    upper = np.asarray(fine['fits']['upper_0.0001_0.01']['phase'])
    lower = np.asarray(fine['fits']['lower_0.0001_0.01']['phase'])
    assert np.max(np.abs(upper - lower)) < .005

    refined = read_run('fem_refined.json')
    shorter = read_run('fem_length8.json')
    initial = read_run('fem_initial.json')
    refinement_change = phase - np.asarray(refined['fits'][reference_key]['phase'])
    domain_change = (np.asarray(refined['fits'][reference_key]['phase']) -
                     np.asarray(shorter['fits'][reference_key]['phase']))
    assert np.max(np.abs(refinement_change)) < .03
    assert np.max(np.abs(domain_change)) < .02

    fit_window_change = (np.asarray(fine['fits']['average_1e-05_0.001']['phase']) - phase)
    assert np.max(np.abs(fit_window_change)) < .005

    # Estimate the oscillation from the computed stresses without imposing
    # epsilon. Terms proportional to r and r^2 account for regular corrections.
    z = np.load('fem_fine.npz')
    r, t = z['r'], (z['upper'] + z['lower']) / 2
    selected = (r > 1e-4) & (r < 1e-2)
    rr = r[selected]
    design = np.column_stack((np.ones_like(rr), np.log(rr), rr, rr**2))
    fitted_epsilon = []
    for case in range(2):
        unwrapped_phase = np.unwrap(np.angle(t[selected, case]))
        coefficients = np.linalg.lstsq(design, unwrapped_phase, rcond=None)[0]
        fitted_epsilon.append(float(coefficients[1]))
    assert np.max(np.abs(np.array(fitted_epsilon) - EPSILON)) < 2e-4

    # Check linear superposition against the mixed term of the energy.
    combined_energy = float(G_COEFFICIENT * abs(kh.sum())**2)
    combined_energy_error = (combined_energy / analytic_energies[2] - 1) * 100
    assert abs(combined_energy_error) < .02

    checks = dict(
        beta=BETA, epsilon=EPSILON, G_over_amplitude_squared=G_COEFFICIENT,
        energies_from_stress_integration=integrated_energies,
        exact_amplitudes=amplitudes.tolist(),
        relative_phase_from_energy_deg=relative_phase_exact,
        relative_phase_from_fem_deg=relative_phase_fem,
        relative_phase_error_deg=relative_phase_error,
        amplitude_error_percent=amplitude_error.tolist(),
        upper_minus_lower_phase_deg=(upper - lower).tolist(),
        refinement_phase_change_deg=refinement_change.tolist(),
        length12_minus_length8_phase_deg=domain_change.tolist(),
        fit_window_phase_change_deg=fit_window_change.tolist(),
        epsilon_from_fem=fitted_epsilon,
        combined_energy_error_percent=combined_energy_error,
        unrounded_final_phases_deg=phase.tolist(),
        unrounded_combined_phase_deg=both_phase,
        mesh_convergence=[dict(
            file=name, nodes=run['nodes'], elements=run['elements'],
            length=run['parameters']['length'],
            maximum_angle_deg=run['parameters']['angle'],
            radial_ratio=run['parameters']['ratio'],
            minimum_radius=run['parameters']['rmin'],
            phase_deg=run['fits'][reference_key]['phase'],
            both_phase_deg=run['fits'][reference_key]['both_phase'],
            amplitudes=run['fits'][reference_key]['amplitude'],
        ) for name, run in [('fem_initial.json', initial), ('fem_refined.json', refined),
                            ('fem_length8.json', shorter), ('fem_fine.json', fine)]],
        all_checks_passed=True,
    )
    Path('checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    Path('results.json').write_text(json.dumps(results, indent=2) + '\n')
    assert json.loads(Path('results.json').read_text()) == results
    assert all(Path(name).is_file() for name in
               ['method.md', 'solve.py', 'postprocess.py', 'checks.json'])
    print(json.dumps(results, indent=2))
    print('All checks passed; details in checks.json.')


if __name__ == '__main__':
    main()
