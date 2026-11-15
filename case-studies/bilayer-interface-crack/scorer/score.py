"""Score one submission of the bilayer interface-crack benchmark against key/reference.json.

Rules stated in the task: eps, the energy release rates and |K| pass within 1 % of the reference value; each phase
angle (psi_deg of the two load cases, psi_both_deg) passes within 0.5 degree (difference taken on the circle). 8 checks. A missing or non-numeric entry fails;
keys not asked for are ignored. Writes score.json and score.md into OUT_DIR. Both hold the reference values: organizer
only. Points by piece of physics are given by rubric.py, from score.json.
Run: python3 -I score.py SUBMISSION_DIR OUT_DIR [VARIANT]   (standard library)
VARIANT: none (default) or no_definitions, the task without the interface-crack definitions
(task/VARIANT_NO_DEFINITIONS.md; reference key/reference_no_definitions.json): G and the amplitudes within 1 %, the
three phase angles at r = h within 0.5 degree, the change of the phase per decade within 0.2 degree. 8 checks.
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF = json.load(open(os.path.join(HERE, '..', 'key', 'reference.json')))
ITEMS = [('eps', None, 'eps', 'rel'),
         ('case A: G E2 h/P^2', 'caseA', 'G_E2_h_over_P2', 'rel'),
         ('case A: |K| sqrt(h)/P', 'caseA', 'Kabs_sqrt_h_over_P', 'rel'),
         ('case A: psi (deg)', 'caseA', 'psi_deg', 'deg'),
         ('case B: G E2 h^3/M^2', 'caseB', 'G_E2_h3_over_M2', 'rel'),
         ('case B: |K| h^(3/2)/M', 'caseB', 'Kabs_h32_over_M', 'rel'),
         ('case B: psi (deg)', 'caseB', 'psi_deg', 'deg'),
         ('both loads, M = P h: psi (deg)', None, 'psi_both_deg', 'deg')]
ITEMS_NO_DEFINITIONS = [('phase change per decade of r, case A (deg)', None, 'phase_change_per_decade_deg', 'abs0.2'),
                        ('case A: G E2 h/P^2', 'caseA', 'G_E2_h_over_P2', 'rel'),
                        ('case A: amplitude sqrt(h)/P', 'caseA', 'amplitude_sqrt_h_over_P', 'rel'),
                        ('case A: phase at r = h (deg)', 'caseA', 'phase_at_h_deg', 'deg'),
                        ('case B: G E2 h^3/M^2', 'caseB', 'G_E2_h3_over_M2', 'rel'),
                        ('case B: amplitude h^(3/2)/M', 'caseB', 'amplitude_h32_over_M', 'rel'),
                        ('case B: phase at r = h (deg)', 'caseB', 'phase_at_h_deg', 'deg'),
                        ('both loads, M = P h: phase at r = h (deg)', None, 'phase_at_h_both_deg', 'deg')]
VARIANTS = {'none': (REF, ITEMS),
            'no_definitions': (json.load(open(os.path.join(HERE, '..', 'key', 'reference_no_definitions.json'))), ITEMS_NO_DEFINITIONS)}


def num(x):
    """A finite JSON number. Booleans, NaN, Infinity and integers too large for a float are not numbers here."""
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return False
    try:
        return math.isfinite(float(x))
    except OverflowError:      # e.g. an integer with 400 digits: valid JSON, not a usable answer
        return False


def get(d, part, key):
    try:
        return d[part][key] if part else d[key]
    except (KeyError, TypeError):
        return None


def score(sub_dir, out_dir, variant='none'):
    ref_all, items = VARIANTS[variant]
    try:
        sub = json.load(open(os.path.join(sub_dir, 'results.json')))
        if not isinstance(sub, dict):
            raise ValueError('results.json is not a JSON object')
        err = None
    except Exception as e:
        sub, err = {}, str(e)[:200]
    rows, passed = [], 0
    for title, part, key, rule in items:
        ref, val = get(ref_all, part, key), get(sub, part, key)
        if not num(val):
            ok = False
        elif rule == 'rel':
            ok = abs(val - ref) <= 0.01 * abs(ref)
        elif rule == 'abs0.2':
            ok = abs(val - ref) <= 0.2
        else:
            ok = abs((val - ref + 180.0) % 360.0 - 180.0) <= 0.5
        rows.append({'item': title, 'key': '%s.%s' % (part, key) if part else key, 'reference': ref,
                     'submitted': val if num(val) else repr(val), 'pass': bool(ok)})
        passed += ok
    os.makedirs(out_dir, exist_ok=True)
    res = {'passed': passed, 'checks': len(items), 'rows': rows, 'read_error': err, 'variant': variant}
    json.dump(res, open(os.path.join(out_dir, 'score.json'), 'w'), indent=1)
    md = ['# Score: %d of %d checks%s' % (passed, len(items), '' if variant == 'none' else ' (variant: %s)' % variant), '']
    if err:
        md += ['results.json could not be read: %s' % err, '']
    md += ['| Item | Reference | Submitted | Pass |', '|---|---|---|---|']
    md += ['| %s | %.6g | %s | %s |' % (r['item'], r['reference'], ('%.6g' % r['submitted']) if num(r['submitted']) else r['submitted'],
                                       'yes' if r['pass'] else '**no**') for r in rows]
    open(os.path.join(out_dir, 'score.md'), 'w').write('\n'.join(md) + '\n')
    return res


if __name__ == '__main__':
    r = score(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else 'none')
    print('%d of %d checks' % (r['passed'], r['checks']), ' | '.join('%s %s' % (x['item'], 'ok' if x['pass'] else 'FAIL') for x in r['rows']))
