"""Controls of the scorer and rubric for the variant without the interface-crack definitions
(task/VARIANT_NO_DEFINITIONS.md, key/reference_no_definitions.json). Right answers must pass; answers with one known
physical error must lose the pieces that the error touches; the tolerance of the new number must sit where the task says.

The wrong answers are the ones this variant is meant to catch, now that the task no longer names the oscillation:
- one_material_tip: the tip treated as that of one material: no change of the phase with r (0), the phase of the
  traction read near the tip (here at r = h/1000, which adds eps ln(1/1000) = +25.5 degrees), the magnitude from G with
  the relation of one material (2.0 % low);
- phase_read_near_tip: the oscillation known, but the phase reported where it was read (r = h/1000), not at r = h;
- reference_length_H: the phase reported at r = H instead of r = h (-2.56 degrees);
- per_decade_sign: the change per decade with the wrong sign;
- per_e_fold: the change per factor e (eps in degrees), not per factor 10;
- eps_itself: the bimaterial constant written in the place of the change per decade;
- one_material_amplitude: the magnitude from G with the relation of one material, all else right.
Run: python3 -I check_controls_no_definitions.py   (standard library; exit 1 if a control misbehaves)
"""
import json, math, os, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score, rubric   # noqa: E402

REF = json.load(open(os.path.join(HERE, '..', 'key', 'reference_no_definitions.json')))
BASE = json.load(open(os.path.join(HERE, '..', 'key', 'reference.json')))
EPS = BASE['eps']
P_TURN, P_G, P_AMP, P_A, P_B, P_BOTH = [p[0] for p in rubric.PIECES_NO_DEFINITIONS]
bad = 0


def report(ok, text):
    global bad
    bad += not ok
    print(('ok   ' if ok else 'FAIL ') + text)


def copy(x):
    return json.loads(json.dumps(x))


def key():
    """A fresh copy each time: the wrong answers below change it in place."""
    return copy({k: REF[k] for k in ('caseA', 'caseB', 'phase_at_h_both_deg', 'phase_change_per_decade_deg')})


def run(obj):
    d = tempfile.mkdtemp()
    try:
        json.dump(obj, open(os.path.join(d, 'results.json'), 'w'))
        s = score.score(d, os.path.join(d, 'out'), 'no_definitions')
        return s, rubric.points(s)
    finally:
        shutil.rmtree(d)


def shift_phases(sub, deg):
    sub['caseA']['phase_at_h_deg'] += deg
    sub['caseB']['phase_at_h_deg'] += deg
    sub['phase_at_h_both_deg'] += deg
    return sub


def scale_amplitudes(sub, f):
    sub['caseA']['amplitude_sqrt_h_over_P'] *= f
    sub['caseB']['amplitude_h32_over_M'] *= f
    return sub


one_mat = math.sqrt(1.0 / math.cosh(math.pi * EPS) ** 2)        # |K| from G with the relation of one material: 0.9798
near = math.degrees(EPS * math.log(1e-3))                         # phase read at r = h/1000 instead of r = h: +25.5
at_H = math.degrees(EPS * math.log(BASE['materials']['h_over_H'] ** -1))   # at r = H = 2h: -2.56

# 1. right answers
s, r = run(key())
report(s['passed'] == 8 and r['points'] == 100, 'key: %d/8, %.0f points' % (s['passed'], r['points']))
sub = key()
for part in ('caseA', 'caseB'):
    sub[part] = {k: round(v, 3) for k, v in sub[part].items()}
sub['phase_at_h_both_deg'] = round(sub['phase_at_h_both_deg'], 2)
sub['phase_change_per_decade_deg'] = round(sub['phase_change_per_decade_deg'], 2)
s, r = run(sub)
report(s['passed'] == 8 and r['points'] == 100, 'key rounded (3 decimals, angles 2): %d/8, %.0f points' % (s['passed'], r['points']))
runs_dir = os.path.join(HERE, '..', 'runs')
for name in ('codex_1', 'kimi_1', 'deepseek_1', 'glm_1'):       # the first round's phase angles, renamed: all inside 0.5 degree
    p = os.path.join(runs_dir, name, 'submission', 'results.json')
    if not os.path.isfile(p):
        continue
    first = json.load(open(p))
    sub = key()
    sub['caseA']['phase_at_h_deg'], sub['caseB']['phase_at_h_deg'] = first['caseA']['psi_deg'], first['caseB']['psi_deg']
    sub['phase_at_h_both_deg'] = first['psi_both_deg']
    s, r = run(sub)
    report(s['passed'] == 8 and r['points'] == 100, 'phase angles of the first-round run %s: %d/8, %.0f points' % (name, s['passed'], r['points']))

# 2. one physical error each: (submission, pieces lost, pieces kept)
wrong = {}
sub = shift_phases(scale_amplitudes(key(), one_mat), near)
sub['phase_change_per_decade_deg'] = 0.0
wrong['one_material_tip'] = (sub, [P_TURN, P_AMP, P_A, P_B, P_BOTH], [P_G])
wrong['phase_read_near_tip'] = (shift_phases(key(), near), [P_A, P_B, P_BOTH], [P_TURN, P_G, P_AMP])
wrong['reference_length_H'] = (shift_phases(key(), at_H), [P_A, P_B, P_BOTH], [P_TURN, P_G, P_AMP])
sub = key(); sub['phase_change_per_decade_deg'] *= -1
wrong['per_decade_sign'] = (sub, [P_TURN], [P_G, P_AMP, P_A, P_B, P_BOTH])
sub = key(); sub['phase_change_per_decade_deg'] = math.degrees(EPS)
wrong['per_e_fold'] = (sub, [P_TURN], [P_G, P_AMP, P_A, P_B, P_BOTH])
sub = key(); sub['phase_change_per_decade_deg'] = EPS
wrong['eps_itself'] = (sub, [P_TURN], [P_G, P_AMP, P_A, P_B, P_BOTH])
wrong['one_material_amplitude'] = (scale_amplitudes(key(), one_mat), [P_AMP], [P_TURN, P_G, P_A, P_B, P_BOTH])
for name, (sub, lost, kept) in wrong.items():
    s, r = run(sub)
    earned = {x['piece']: x['earned'] for x in r['rows']}
    ok = all(earned[p] == 0 for p in lost) and all(earned[p] > 0 for p in kept)
    report(ok, '%s: %d/8, %.0f points; lost: %s' % (name, s['passed'], r['points'], '; '.join(p[:48] for p in r['missing'])))

# 3. the tolerance of the change per decade: 0.2 degree
for d, want in ((0.19, 8), (0.21, 7), (-0.19, 8), (-0.21, 7)):
    sub = key(); sub['phase_change_per_decade_deg'] += d
    s, _ = run(sub)
    report(s['passed'] == want, 'edge: change per decade %+.2f degree gives %d/8 (expected %d)' % (d, s['passed'], want))

# 4. a right answer in the form of the first round (the definitions' keys): only the two energy release rates, whose
#    names the two versions share, are read; every other entry is missing here and fails
first_form = {'eps': EPS, 'caseA': {'G_E2_h_over_P2': BASE['caseA']['G_E2_h_over_P2'], 'Kabs_sqrt_h_over_P': BASE['caseA']['Kabs_sqrt_h_over_P'],
                                   'psi_deg': BASE['caseA']['psi_deg']},
              'caseB': {'G_E2_h3_over_M2': BASE['caseB']['G_E2_h3_over_M2'], 'Kabs_h32_over_M': BASE['caseB']['Kabs_h32_over_M'],
                        'psi_deg': BASE['caseB']['psi_deg']}, 'psi_both_deg': BASE['psi_both_deg']}
s, r = run(first_form)
report(s['passed'] == 2 and r['points'] == 20, 'a right answer under the first round\'s keys: %d/8, %.0f points (the two G only)' % (s['passed'], r['points']))

print('note: one material near the tip: amplitude factor %.4f, phase read at h/1000 shifts by %+.2f degrees; at r = H by %+.2f'
      % (one_mat, near, at_H))
print('ALL CONTROLS BEHAVE' if bad == 0 else '%d CONTROL(S) MISBEHAVE' % bad)
sys.exit(1 if bad else 0)
