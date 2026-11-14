"""Physics rubric of the trapped-front benchmark: points by piece of physics, on top of score.py.

Teng, 13 Nov 2026, after the first four runs: "answer close is one thing. missing a physics is another. If one agent
misses a physics, they need to lose points"; "We just need to update the rubric".

score.py still decides whether a number is right, by the rule stated in the task (within 1% of the reference or within
0.005). This file weighs the right numbers by the physics they test:

- The response, 70 points: seven pieces of physics, 10 points each. Each is tested by one or two numbers of Part A,
  which isolate it, and is earned only if all of them pass. A run that misses a piece loses that piece's points.
- The outcome, 30 points: the equilibrium front and the stress intensity factors along it (Part B). It needs every
  piece of the response, so it is counted only in a run that has all seven; then each part earns its points in
  proportion to its numbers that pass (a slip in one number costs little, a missing piece of physics costs all).

A number that is close but outside the tolerance earns nothing, here as in the count. The load level is a Part B number
but tests the energy relation alone, so it sits with the response.

The points (10 a piece, 30 for the outcome) and the rule that the outcome counts only with a complete response were
proposed by us and approved by Teng on 13 Nov 2026 ("please document. I think it is ok."). The runs of 13 Nov 2026 were
told only the pass rule; the rubric was written after them. The rubric in words: RUBRIC.md. Organizer only.

usage: python3 -I rubric.py SCORE_DIR       reads SCORE_DIR/score.json, writes rubric.json and rubric.md there
"""
import json, os, sys

# (piece of physics, points, the numbers that test it)
RESPONSE = [
    ('Tension: a part of the front that lags is loaded more', 10, ['C1', 'S1']),
    ('Energy release rate in mixed mode and the load level of equilibrium', 10, ['K10_scaled']),
    ("Shear across the front: its own weight, which depends on Poisson's ratio", 10, ['C2']),
    ("Shear along the front: its own weight, which depends on Poisson's ratio", 10, ['C3']),
    ('The slope of the front turns shear along it into shear across it', 10, ['S2']),
    ('The slope of the front turns shear across it into shear along it', 10, ['S3']),
    ('The energy release rate of a wavy front, with no term in the slope', 10, ['CG', 'SG']),
]
OUTCOME = [
    ('The equilibrium front in the tough bands', 10, ['front_peak_to_peak'] + ['f[%d]' % k for k in range(16)]),
    ('The stress intensity factors along that front', 20,
     ['%s[%d]' % (t, k) for t in ('kappa1', 'kappa2', 'kappa3') for k in range(16)] + ['kappa3_max', 's_kappa3_max']),
]
_all = [i for _, _, items in RESPONSE + OUTCOME for i in items]
assert len(_all) == 76 and len(set(_all)) == 76 and sum(p for _, p, _ in RESPONSE + OUTCOME) == 100


def points(score):
    """score: the dict score.py writes (score.json). Returns the rubric result as a dict."""
    ok = {r['item']: bool(r['pass']) for r in score['rows']}
    assert set(ok) == set(_all), 'score.json does not hold the 76 checks of this benchmark'
    rows = []
    for title, pts, items in RESPONSE:
        n = sum(ok[i] for i in items)
        rows.append({'part': 'response', 'piece': title, 'points': pts, 'numbers': len(items), 'passed': n,
                     'earned': float(pts) if n == len(items) else 0.0})
    whole = all(r['earned'] == r['points'] for r in rows)
    for title, pts, items in OUTCOME:
        n = sum(ok[i] for i in items)
        rows.append({'part': 'outcome', 'piece': title, 'points': pts, 'numbers': len(items), 'passed': n,
                     'earned': pts * n / len(items) if whole else 0.0})
    return {'points': sum(r['earned'] for r in rows), 'of': 100, 'response_complete': whole,
            'missing': [r['piece'] for r in rows if r['part'] == 'response' and r['earned'] == 0], 'rows': rows,
            'checks_passed': score['passed'], 'checks': score['checks']}


def markdown(res):
    md = ['# Physics rubric: %.0f of 100 points (%d of %d numbers within the tolerance)' % (res['points'], res['checks_passed'], res['checks']), '',
          '| Piece of physics | Numbers passed | Points |', '|---|---|---|']
    for r in res['rows']:
        md.append('| %s | %d / %d | %s / %d |' % (r['piece'], r['passed'], r['numbers'], ('%.1f' % r['earned']).rstrip('0').rstrip('.'), r['points']))
    if not res['response_complete']:
        md += ['', 'The outcome is not counted: it needs every piece of the response, and %d %s missing.' %
               (len(res['missing']), 'is' if len(res['missing']) == 1 else 'are')]
    return '\n'.join(md) + '\n'


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    d = sys.argv[1]
    res = points(json.load(open(os.path.join(d, 'score.json'))))
    json.dump(res, open(os.path.join(d, 'rubric.json'), 'w'), indent=1)
    open(os.path.join(d, 'rubric.md'), 'w').write(markdown(res))
    print('%.0f of 100 points; response complete: %s; missing: %d pieces' % (res['points'], res['response_complete'], len(res['missing'])))
