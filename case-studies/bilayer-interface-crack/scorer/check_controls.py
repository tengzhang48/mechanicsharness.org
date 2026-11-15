"""Controls of the bilayer scorer and rubric: right answers must pass; answers with a known physical error must lose
the checks and the pieces of the rubric that the error touches; broken files must fail without stopping the scorer.

Each wrong answer is built from the key's closed forms and phase angles with one change, the kind of change a solver
could make: the phase angle of one material used for two; the paper's printed table entry used for the phase angle; the
lower layer's thickness used as the reference length; the one-material relation between G and |K|; plane stress; the sign
convention of K2; the energy of the upper arm only; the two load cases exchanged; the two loadings combined by averaging
their phase angles.

Run: python3 -I check_controls.py      (standard library; reads ../key/reference.json; exit 1 if a control misbehaves)
"""
import cmath, json, math, os, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score   # noqa: E402
import rubric  # noqa: E402

REF = json.load(open(os.path.join(HERE, '..', 'key', 'reference.json')))
M = REF['materials']
E1, NU1, NU2, ETA = M['E1_over_E2'], M['nu1'], M['nu2'], M['h_over_H']      # E2 = 1, h = 1, H = 1/ETA
CHI = REF['M_over_Ph_both']
P_EPS, P_G, P_K, P_A, P_B, P_BOTH = [p[0] for p in rubric.PIECES]
bad = 0


def report(ok, text):
    global bad
    bad += not ok
    print(('ok   ' if ok else 'FAIL ') + text)


def constants(plane_stress=False):
    """Dundurs' alpha, the bimaterial constant eps, c1 + c2 and the two effective moduli, for plane strain or plane stress."""
    mu1, mu2 = E1 / (2 * (1 + NU1)), 1.0 / (2 * (1 + NU2))
    kap = (lambda nu: (3 - nu) / (1 + nu)) if plane_stress else (lambda nu: 3 - 4 * nu)
    k1, k2 = kap(NU1), kap(NU2)
    eps = math.log((k1 / mu1 + 1 / mu2) / (k2 / mu2 + 1 / mu1)) / (2 * math.pi)
    Eb1, Eb2 = (E1, 1.0) if plane_stress else (E1 / (1 - NU1 ** 2), 1.0 / (1 - NU2 ** 2))
    return (Eb1 - Eb2) / (Eb1 + Eb2), eps, (k1 + 1) / mu1 + (k2 + 1) / mu2, Eb1, Eb2


def gamma_deg(alpha):
    sig = (1 + alpha) / (1 - alpha)
    A = 1 / (1 + sig * (4 * ETA + 6 * ETA ** 2 + 3 * ETA ** 3))
    I = 1 / (12 * (1 + sig * ETA ** 3))
    return math.degrees(math.asin(6 * sig * ETA ** 2 * (1 + ETA) * math.sqrt(A * I)))


def psi_both_of(sub):
    """Phase angle under both loads, M = CHI P h, from a submission's two magnitudes and phase angles."""
    kA = sub['caseA']['Kabs_sqrt_h_over_P'] * cmath.exp(1j * math.radians(sub['caseA']['psi_deg']))
    kB = sub['caseB']['Kabs_h32_over_M'] * cmath.exp(1j * math.radians(sub['caseB']['psi_deg']))
    return math.degrees(cmath.phase(kA + CHI * kB))


def key():
    return {'eps': REF['eps'], 'caseA': {k: REF['caseA'][k] for k in ('G_E2_h_over_P2', 'Kabs_sqrt_h_over_P', 'psi_deg')},
            'caseB': {k: REF['caseB'][k] for k in ('G_E2_h3_over_M2', 'Kabs_h32_over_M', 'psi_deg')}, 'psi_both_deg': REF['psi_both_deg']}


def copy(sub):
    return json.loads(json.dumps(sub))


def with_omega(omega):
    """The key with the phase angle omega of the source's solution replaced: psi_A = omega, psi_B = omega + gamma - 90."""
    s = key()
    s['caseA']['psi_deg'] = omega
    s['caseB']['psi_deg'] = omega + gamma_deg(REF['alpha']) - 90
    s['psi_both_deg'] = psi_both_of(s)
    return s


def rounded(sub, nd):
    r = lambda v: round(v, nd)
    return {k: ({kk: r(vv) for kk, vv in v.items()} if isinstance(v, dict) else r(v)) for k, v in sub.items()}


tmp = tempfile.mkdtemp()


def run(name, obj):
    d = os.path.join(tmp, name)
    os.makedirs(d)
    if obj is not None:
        with open(os.path.join(d, 'results.json'), 'w') as fh:
            fh.write(obj if isinstance(obj, str) else json.dumps(obj))
    s = score.score(d, d + '_out')
    return s, rubric.points(s)


try:
    good = key()
    alpha, eps, csum, Eb1, Eb2 = constants()
    # 0. the key is consistent with its own closed forms
    report(abs(eps - REF['eps']) < 1e-12 and abs(alpha - REF['alpha']) < 1e-12, 'eps and alpha recomputed here equal the key: %.6f, %.4f' % (eps, alpha))
    report(abs(REF['caseB']['psi_deg'] - (REF['caseA']['psi_deg'] + gamma_deg(alpha) - 90)) < 0.001,
           'the two phase angles of the key differ by gamma - 90 degrees (exact relation, from energy) within 0.001 degree: %.6f'
           % (REF['caseB']['psi_deg'] - (REF['caseA']['psi_deg'] + gamma_deg(alpha) - 90)))
    report(abs(psi_both_of(good) - REF['psi_both_deg']) < 1e-9, 'the phase angle under both loads follows from the two complex K of the key')
    for case, gk, kk in (('caseA', 'G_E2_h_over_P2', 'Kabs_sqrt_h_over_P'), ('caseB', 'G_E2_h3_over_M2', 'Kabs_h32_over_M')):
        k_from_g = math.sqrt(16 * math.cosh(math.pi * eps) ** 2 * good[case][gk] / csum)
        report(abs(k_from_g / good[case][kk] - 1) < 1e-12, '%s: |K| of the key follows from its G by the two-material relation' % case)
    if 'first_route' in REF:
        d = max(abs(REF['first_route'][k]) for k in ('diff_psiA_deg', 'diff_psiB_deg', 'diff_psi_both_deg'))
        report(d < 0.05, 'the first route\'s three phase angles are within %.3f degree of the key\'s, which are the second route\'s (0.05 allowed; the tolerance of the task is 0.5)' % d)
    else:
        report(False, 'the key holds no check by the first route')

    # 1. right answers
    for name, sub in (('key', good), ('key_3_decimals', rounded(good, 3)), ('psi_B_plus_360', None)):
        if sub is None:
            sub = copy(good)
            sub['caseB']['psi_deg'] += 360.0
        s, r = run(name, sub)
        report(s['passed'] == 8 and r['points'] == 100, '%s: %d/8, %.0f points' % (name, s['passed'], r['points']))
    if 'first_route' in REF:
        sub = copy(good)
        sub['caseA']['psi_deg'], sub['caseB']['psi_deg'] = REF['first_route']['psiA_deg'], REF['first_route']['psiB_deg']
        sub['psi_both_deg'] = REF['first_route']['psi_both_deg']
        s, r = run('first_route_phase_angles', sub)
        report(s['passed'] == 8 and r['points'] == 100, 'the first route\'s phase angles as the answer: %d/8, %.0f points' % (s['passed'], r['points']))

    # 2. wrong physics: name -> (submission, pieces that must earn nothing, pieces that must be earned)
    wrong = {}
    wrong['one_material_phase_angle'] = (with_omega(52.1 - 3 * ETA), [P_A, P_B, P_BOTH], [P_EPS, P_G, P_K])
    wrong['printed_table_entry'] = (with_omega(REF['record_printed_table']['omega_deg']), [P_A, P_B, P_BOTH], [P_EPS, P_G, P_K])
    s3 = copy(good)                                          # the lower layer's thickness H as the reference length
    shift = math.degrees(eps * math.log(1 / ETA))
    for k in ('caseA', 'caseB'):
        s3[k]['psi_deg'] += shift
    s3['psi_both_deg'] += shift
    wrong['reference_length_H'] = (s3, [P_A, P_B, P_BOTH], [P_EPS, P_G, P_K])
    s4 = copy(good)                                          # |K| from G without the factor for two materials
    s4['caseA']['Kabs_sqrt_h_over_P'] /= math.cosh(math.pi * eps)
    s4['caseB']['Kabs_h32_over_M'] /= math.cosh(math.pi * eps)
    wrong['one_material_K_from_G'] = (s4, [P_K], [P_EPS, P_G, P_A, P_B, P_BOTH])
    a_ps, eps_ps, csum_ps, Eb1_ps, Eb2_ps = constants(plane_stress=True)
    assert abs(NU1 - NU2) < 1e-12, 'the plane-stress control below assumes equal Poisson ratios (G scales by 1/(1 - nu^2))'
    s5 = copy(good)                                          # plane stress constants in eps, G and |K|
    s5['eps'] = eps_ps
    for case, gk, kk in (('caseA', 'G_E2_h_over_P2', 'Kabs_sqrt_h_over_P'), ('caseB', 'G_E2_h3_over_M2', 'Kabs_h32_over_M')):
        s5[case][gk] = good[case][gk] / (1 - NU2 ** 2)
        s5[case][kk] = math.sqrt(16 * math.cosh(math.pi * eps_ps) ** 2 * s5[case][gk] / csum_ps)
    wrong['plane_stress'] = (s5, [P_EPS, P_G], [P_A, P_B, P_BOTH])
    s6 = copy(good)                                          # the other sign convention for K2
    for k in ('caseA', 'caseB'):
        s6[k]['psi_deg'] *= -1
    s6['psi_both_deg'] *= -1
    wrong['sign_of_K2'] = (s6, [P_A, P_B, P_BOTH], [P_EPS, P_G, P_K])
    s7 = copy(good)                                          # energy of the upper arm only
    s7['caseA']['G_E2_h_over_P2'] = 1 / (2 * Eb1)
    s7['caseB']['G_E2_h3_over_M2'] = 6 / Eb1
    for case, gk, kk in (('caseA', 'G_E2_h_over_P2', 'Kabs_sqrt_h_over_P'), ('caseB', 'G_E2_h3_over_M2', 'Kabs_h32_over_M')):
        s7[case][kk] = math.sqrt(16 * math.cosh(math.pi * eps) ** 2 * s7[case][gk] / csum)
    s7['psi_both_deg'] = psi_both_of(s7)
    wrong['upper_arm_energy_only'] = (s7, [P_G, P_K], [P_EPS, P_A, P_B])
    s8 = copy(good)                                          # the two load cases exchanged
    s8['caseA'] = {'G_E2_h_over_P2': good['caseB']['G_E2_h3_over_M2'], 'Kabs_sqrt_h_over_P': good['caseB']['Kabs_h32_over_M'], 'psi_deg': good['caseB']['psi_deg']}
    s8['caseB'] = {'G_E2_h3_over_M2': good['caseA']['G_E2_h_over_P2'], 'Kabs_h32_over_M': good['caseA']['Kabs_sqrt_h_over_P'], 'psi_deg': good['caseA']['psi_deg']}
    wrong['cases_exchanged'] = (s8, [P_G, P_K, P_A, P_B], [P_EPS])
    s9 = copy(good)                                          # the two loadings combined by averaging the phase angles
    s9['psi_both_deg'] = 0.5 * (good['caseA']['psi_deg'] + good['caseB']['psi_deg'])
    wrong['phase_angles_averaged'] = (s9, [P_BOTH], [P_EPS, P_G, P_K, P_A, P_B])
    for name, (sub, lost, kept) in wrong.items():
        s, r = run(name, sub)
        earned = {x['piece']: x['earned'] for x in r['rows']}
        ok = all(earned[p] == 0 for p in lost) and all(earned[p] > 0 for p in kept) and r['points'] < 100
        report(ok, '%s: %d/8, %.0f points; lost: %s' % (name, s['passed'], r['points'], '; '.join(p.split(':')[0][:44] for p in r['missing'])))
    d_tab = REF['caseA']['psi_deg'] - REF['record_printed_table']['omega_deg']
    print('note: the printed table entry is %.2f degrees from the key\'s phase angle under the force (the check allows 0.5)' % d_tab)

    # 3. the tolerance edge: 1 % on a magnitude, 0.5 degree on a phase angle
    for part, k, f_ok, f_bad in (('caseA', 'G_E2_h_over_P2', 1.009, 1.011), ('caseB', 'Kabs_h32_over_M', 0.991, 0.989)):
        for factor, want in ((f_ok, 8), (f_bad, 7)):
            sub = copy(good)
            sub[part][k] *= factor
            s, _ = run('edge_%s_%s' % (k, factor), sub)
            report(s['passed'] == want, 'edge: %s times %.3f gives %d/8 (expected %d)' % (k, factor, s['passed'], want))
    for d, want in ((0.45, 8), (0.55, 7), (-0.45, 8), (-0.55, 7)):
        sub = copy(good)
        sub['caseA']['psi_deg'] += d
        s, _ = run('edge_psi_%s' % d, sub)
        report(s['passed'] == want, 'edge: psi of case A shifted by %+.2f degrees gives %d/8 (expected %d)' % (d, s['passed'], want))

    # 4. broken files
    huge = json.dumps(good).replace(json.dumps(good['eps']), '1' + '0' * 400, 1)
    nan = json.dumps(good).replace(json.dumps(good['caseA']['psi_deg']), 'NaN', 1)
    assert huge != json.dumps(good) and nan.count('NaN') == 1
    strings = copy(good)
    strings['caseB'] = {k: str(v) for k, v in good['caseB'].items()}
    booleans = copy(good)
    booleans['eps'] = True
    no_a = copy(good)
    del no_a['caseA']
    extra = copy(good)
    extra['note'] = 'extra key'
    extra['caseA']['omega'] = 50.0
    for name, obj, want in (('empty_folder', None, 0), ('malformed_json', '{"eps": ', 0), ('not_object', '[1, 2]', 0),
                            ('huge_integer_eps', huge, 7), ('nan_psi', nan, 7), ('strings_in_caseB', strings, 5),
                            ('boolean_eps', booleans, 7), ('caseA_missing', no_a, 5), ('extra_keys_ignored', extra, 8)):
        try:
            s, r = run(name, obj)
            report(s['passed'] == want, '%s: %d/8 (expected %d), %.0f points' % (name, s['passed'], want, r['points']))
        except Exception as e:      # a crash of the scorer is a failed control
            report(False, '%s: the scorer crashed: %s: %s' % (name, type(e).__name__, str(e)[:80]))
finally:
    shutil.rmtree(tmp)
print('FAILED %d' % bad if bad else 'ALL CONTROLS BEHAVE')
sys.exit(1 if bad else 0)
