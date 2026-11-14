"""Score one submission of the trapped-front benchmark against key/reference.json.

A number passes if |x - x_ref| <= max(0.01 |x_ref|, 0.005) (the rule stated in the task). Checks: Part A 8 numbers;
Part B K10_scaled, front_peak_to_peak, 16 values each of f, kappa1, kappa2, kappa3, kappa3_max, s_kappa3_max (68).
A missing or non-numeric entry fails its check; keys not asked for are ignored. Writes score.json and score.md into
OUT_DIR. Both files contain the reference values: organizer only, never shown to an agent during a run.
Run: python3 -I score.py SUBMISSION_DIR OUT_DIR   (standard library only)
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF = json.load(open(os.path.join(HERE, '..', 'key', 'reference.json')))
GROUPS = [('Part A: cosine-front coefficients', 'partA', ['C1', 'S1', 'C2', 'S2', 'C3', 'S3', 'CG', 'SG']),
          ('Part B: loading level and front amplitude', 'partB', ['K10_scaled', 'front_peak_to_peak']),
          ('Part B: front shape f(s)', 'partB', [('f', k) for k in range(16)]),
          ('Part B: kappa1(s)', 'partB', [('kappa1', k) for k in range(16)]),
          ('Part B: kappa2(s)', 'partB', [('kappa2', k) for k in range(16)]),
          ('Part B: kappa3(s)', 'partB', [('kappa3', k) for k in range(16)]),
          ('Part B: maximum of kappa3', 'partB', ['kappa3_max', 's_kappa3_max'])]


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
        v = d[part]
        if isinstance(key, tuple):
            v = v[key[0]]
            if not isinstance(v, list) or len(v) != 16:
                return None
            return v[key[1]]
        return v[key]
    except (KeyError, TypeError, IndexError):
        return None


def score(sub_dir, out_dir):
    path = os.path.join(sub_dir, 'results.json')
    try:
        sub = json.load(open(path))
        if not isinstance(sub, dict):
            raise ValueError('results.json is not a JSON object')
        err = None
    except Exception as e:  # unreadable submission: every check fails
        sub, err = {}, str(e)[:200]
    rows, total, passed = [], 0, 0
    summary = []
    for title, part, keys in GROUPS:
        gp = 0
        for key in keys:
            ref = get(REF, part, key)
            val = get(sub, part, key)
            ok = num(val) and abs(val - ref) <= max(0.01 * abs(ref), 0.005)
            name = '%s[%d]' % key if isinstance(key, tuple) else key
            rows.append({'group': title, 'item': name, 'reference': ref, 'submitted': val if num(val) else repr(val),
                         'pass': bool(ok)})
            gp += ok
        summary.append({'group': title, 'passed': gp, 'checks': len(keys)})
        total += len(keys)
        passed += gp
    os.makedirs(out_dir, exist_ok=True)
    res = {'passed': passed, 'checks': total, 'groups': summary, 'rows': rows, 'read_error': err}
    json.dump(res, open(os.path.join(out_dir, 'score.json'), 'w'), indent=1)
    md = ['# Score: %d of %d checks' % (passed, total), '']
    if err:
        md += ['results.json could not be read: %s' % err, '']
    md += ['| Group | Passed |', '|---|---|'] + ['| %s | %d / %d |' % (g['group'], g['passed'], g['checks']) for g in summary]
    md += ['', '| Item | Reference | Submitted | Pass |', '|---|---|---|---|']
    md += ['| %s | %.6g | %s | %s |' % (r['item'], r['reference'], ('%.6g' % r['submitted']) if num(r['submitted']) else r['submitted'],
                                       'yes' if r['pass'] else '**no**') for r in rows]
    open(os.path.join(out_dir, 'score.md'), 'w').write('\n'.join(md) + '\n')
    return res


if __name__ == '__main__':
    r = score(sys.argv[1], sys.argv[2])
    print('%d of %d checks' % (r['passed'], r['checks']), ' | '.join('%s %d/%d' % (g['group'], g['passed'], g['checks']) for g in r['groups']))
