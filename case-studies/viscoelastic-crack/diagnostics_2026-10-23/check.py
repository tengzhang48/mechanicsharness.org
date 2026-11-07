"""Focused source-based diagnostics; no solver submission or grade is changed.

Run: python3 check.py
Requires NumPy, SciPy, mpmath. Original submissions are read-only inputs.

Main source: Hui and Long (2026), equations 12, 14, 17, 51, 59.
The causal control uses the exact mixed-boundary stress field, including its
homogeneous Griffith mode, but replaces zero true singular coefficient with
Azure DeepSeek's finite-window regression condition. All other continuum
conditions remain exact. This isolates a selection error from mesh errors.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path('<pilot-folder>')
OUT = Path(__file__).resolve().parent / 'results.json'
PI = np.pi
LAM = 0.1
DELTA = 0.2


def loading(epsilon):
    return 2 / PI * np.arccos(np.exp(-epsilon / 2))


def slope(r, epsilon):
    den = -np.expm1(-epsilon)
    c = np.exp(-r / 2) * np.sqrt(-np.expm1(-(epsilon-r)) / den)
    return 2 / PI * (np.log1p(c) - 0.5 * np.log(-np.expm1(-r) / den))


def memory(epsilon, speed):
    return speed * quad(
        lambda t: np.exp(-t) * slope(speed*t, epsilon),
        0, min(45., epsilon/speed), epsabs=2e-12, epsrel=2e-12)[0]


def opening(epsilon, speed):
    return PI * loading(epsilon)**2 / 2 - (1-LAM) * memory(epsilon, speed)


def fit_coefficient(x, stress, extra_regular_terms=0):
    powers = [-.5, 0., .5] + [1.5+j for j in range(extra_regular_terms)]
    design = np.column_stack([x**power for power in powers])
    return float(np.linalg.lstsq(design, stress, rcond=None)[0][0])


def bounded_traction(x, epsilon):
    return 1 - 2/PI * np.arctan(
        np.sqrt(-np.expm1(-x)) / np.sqrt(np.expm1(epsilon)))


def fit_load_shift(epsilon, panel_nodes, extra_regular_terms=0):
    # Exact same front-panel coordinates and fixed fit window as Azure.
    x = 14 * (1-np.cos(np.linspace(0, PI, panel_nodes+1))) / 2
    x = x[(x > 0.02) & (x < 0.4)]
    bounded = bounded_traction(x, epsilon)
    griffith = 1 / np.sqrt(-np.expm1(-x))
    return -fit_coefficient(x, bounded, extra_regular_terms) / fit_coefficient(
        x, griffith, extra_regular_terms)


def griffith_opening(epsilon, speed):
    # u=epsilon*t^2 regularizes the homogeneous mode's front singularity.
    return quad(
        lambda t: 2*epsilon*t / np.sqrt(np.expm1(epsilon*t*t)) *
        (1-(1-LAM)*np.exp(-epsilon*(1-t*t)/speed)),
        0, 1, epsabs=2e-11, epsrel=2e-11)[0]


def main():
    result = {
        'checked_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Source-based diagnostic calculations, not an independent scientific grade.',
        'source_equations': 'Hui and Long main equations 12, 14, 17, 51, 59',
        'original_files_changed': False,
        'input_sha256': {},
    }
    paths = ['private/<source paper PDF, not included>', 'private/REFERENCE.md']
    names = ['P01_deepseek', 'P04_deepseek_direct_1017', 'P03_gpt61_sol']
    cases = {}
    for name in names:
        for filename in ['ANSWER.md', 'solve.py', 'results.json']:
            paths.append(f'runs/{name}/workspace/{filename}')
        cases[name] = json.loads((ROOT/f'runs/{name}/workspace/results.json').read_text())
    for rel in paths:
        result['input_sha256'][rel] = hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()

    relation_checks = []
    for name, data in cases.items():
        for case in data['cases']:
            if not case['admissible']:
                continue
            actual = case['a']
            necessary = loading(case['epsilon'])
            relation_checks.append({
                'run': name, 'V': case['V'], 'epsilon': case['epsilon'],
                'a': actual, 'a_required_by_bounded_stress': float(necessary),
                'true_front_singular_coefficient': float(actual-necessary),
                'relative_load_mismatch': float((actual-necessary)/necessary),
            })
    result['necessary_load_zone_relation'] = relation_checks

    isolation = []
    for speed, reported, physical in [
        (0.001, .577, .5447055550934302),
        (0.1, .7686, .7407809316699946),
    ]:
        for nodes in [320, 400]:
            epsilon = brentq(
                lambda e: opening(e, speed) + fit_load_shift(e, nodes) *
                griffith_opening(e, speed) - PI*DELTA/2,
                .4, 1.2, xtol=1e-12)
            isolation.append({
                'V': speed, 'front_panel_nodes': nodes,
                'true_epsilon': physical, 'reported_azure_epsilon': reported,
                'epsilon_with_only_fit_condition_changed': float(epsilon),
                'a_with_only_fit_condition_changed': float(
                    loading(epsilon)+fit_load_shift(epsilon, nodes)),
                'fraction_of_reported_zone_error_reproduced': float(
                    (epsilon-physical)/(reported-physical)),
            })
    result['azure_fit_isolation'] = isolation
    corrections = []
    for speed, physical in [(0.001, .5447055550934302),
                            (0.1, .7407809316699946)]:
        for terms in range(5):
            epsilon = brentq(
                lambda e: opening(e, speed) + fit_load_shift(e, 400, terms) *
                griffith_opening(e, speed) - PI*DELTA/2,
                .4, 1.1, xtol=1e-12)
            corrections.append({
                'V': speed, 'extra_regular_terms_in_fit': terms,
                'added_powers': [1.5+j for j in range(terms)],
                'predicted_epsilon': float(epsilon),
                'relative_zone_error': float((epsilon-physical)/physical),
                'unchanged': 'Exact continuum stress and opening, original fit interval, original front-panel coordinates',
            })
    result['azure_omitted_regular_terms_correction_control'] = corrections
    result['why_the_fit_is_biased'] = {
        'original_basis': ['X**(-1/2)', '1', 'X**(1/2)'],
        'true_bounded_front_expansion': 'T2 = 1 - 2/(pi*L)*sqrt(X) + [1/(2*pi*L)+2/(3*pi*L**3)]*X**(3/2) + O(X**(5/2)); L=sqrt(exp(epsilon)-1)',
        'mechanism': 'Omitted regular curvature projects into the fitted X**(-1/2) coefficient. Forcing that regression coefficient to zero changes the physical load and zone length.',
        'mesh_refinement_limit': 'A fixed fit interval retains this bias even with exact continuum fields.',
        'speed_dependence': 'The omitted curvature is larger for short low-speed zones and smaller for the long high-speed zone.',
        'leading_omitted_term_examples': [
            {'epsilon': epsilon,
             'coefficient_of_X_to_three_halves': float(
                 1/(2*PI*np.sqrt(np.expm1(epsilon))) +
                 2/(3*PI*np.expm1(epsilon)**1.5))}
            for epsilon in [.577, .7686, 6.34]],
    }
    result['isolation_limitations'] = [
        'This isolates the regression condition in the continuum problem; it does not replay all numerical discretization errors.',
        'Errors remaining after this control can include the mesh, constitutive averaging, field reconstruction and rounded selected states.',
        'The finite-window coefficient is a regression parameter, not the true asymptotic stress-intensity coefficient.',
    ]

    sensitivity = []
    for scale in [1, .3, .1, .03, .01]:
        x = np.geomspace(.02*scale, .4*scale, 200)
        sensitivity.append({
            'fit_window': [float(x[0]), float(x[-1])],
            'true_singular_coefficient': 0.,
            'fitted_singular_coefficient': fit_coefficient(x, bounded_traction(x, .577)),
        })
    result['fit_on_exact_bounded_stress_window_sensitivity'] = sensitivity

    e0 = -2*np.log(np.cos(PI*np.sqrt(DELTA)/2))
    slow = []
    for speed in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
        epsilon = brentq(lambda e: opening(e, speed)-PI*DELTA/2,
                         e0, 2*e0, xtol=1e-14)
        loss = 2*(1-LAM)*memory(epsilon, speed)/(PI*DELTA)
        leading = 2*(1-LAM)*speed/(PI**2*DELTA) * (
            np.log(4*(-np.expm1(-e0))/speed)+np.euler_gamma)
        slow.append({'V': speed, 'exact_scalar_loss': float(loss),
                     'gpt_leading_expansion': float(leading),
                     'relative_error': float(abs(leading-loss)/loss),
                     'remainder_over_V': float((loss-leading)/speed)})
    result['gpt_optional_slow_speed_expansion'] = slow

    mp.mp.dps = 60
    lam, delta, pi = mp.mpf('.1'), mp.mpf('.2'), mp.pi

    def k(v):
        return 2*v/pi**2 * (mp.digamma(mp.mpf('.5')+1/v)-mp.digamma(mp.mpf('.5')))

    def h_infinite(v):
        return pi/2*(1-(1-lam)*k(v))

    terminal = mp.findroot(lambda v: 1-delta-(1-lam)*k(v), (12, 15))

    def h_finite(e, v):
        a = 2/pi*mp.acos(mp.exp(-e/2))
        den = -mp.expm1(-e)

        def integrand(t):
            if t == 0 or t == 1:
                return mp.mpf(0)
            r = e*t*t
            c = mp.exp(-r/2)*mp.sqrt(-mp.expm1(-(e-r))/den)
            s = 2/pi*(mp.log1p(c)-mp.log(-mp.expm1(-r)/den)/2)
            return 2*e*t*mp.exp(-r/v)*s

        integral = mp.quad(integrand, [0, mp.mpf('.01'), mp.mpf('.1'),
                                      mp.mpf('.5'), mp.mpf('.9'), 1])
        return pi*a*a/2-(1-lam)*integral

    result['gpt_endpoint_checks'] = {
        'Vstar_from_source_condition_60_digit_arithmetic': mp.nstr(terminal, 40),
        'minus_Hinfinite_derivative': mp.nstr(-mp.diff(h_infinite, terminal), 30),
        'gap_divided_by_2_exp_minus_epsilon_over_2': [
            {'epsilon': e, 'ratio_expected_to_approach_one': mp.nstr(
                (h_infinite(terminal)-h_finite(mp.mpf(e), terminal)) /
                (2*mp.exp(-mp.mpf(e)/2)), 25)}
            for e in [20, 40, 60, 80, 100]],
    }
    result['terminal_depends_on_fracture_threshold_at_fixed_lambda'] = [
        {'lambda': '.1', 'delta_bar': str(dd), 'terminal_speed': mp.nstr(
            mp.findroot(lambda v: 1-dd-(1-lam)*k(v), guess), 30)}
        for dd, guess in [(mp.mpf('.15'), 30), (mp.mpf('.2'), 13), (mp.mpf('.3'), 5)]]

    result['gpt_check_scope'] = [
        'No discrepancy found in these targeted source-relation and asymptotic checks.',
        'Numerical asymptotic convergence is evidence, not a uniform error proof.',
        'The prior external grade is unchanged; no new total score is assigned.',
        'These calculations do not establish performance on a new blind problem or identify internal reasoning causes.',
    ]
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(f'Saved {OUT}')
    print('Azure zone-error fraction reproduced:', [
        round(x['fraction_of_reported_zone_error_reproduced'], 4)
        for x in isolation if x['front_panel_nodes'] == 400])
    print('Source terminal speed:', mp.nstr(terminal, 20))
    print('GPT slow-expansion relative error at V=1e-7:', slow[-1]['relative_error'])


if __name__ == '__main__':
    main()
