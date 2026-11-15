"""Physics rubric of the bilayer interface-crack benchmark: points by piece of physics, on top of score.py.

The principle is Teng's, set on the trapped-front benchmark on 13 Nov 2026: "answer close is one thing. missing a
physics is another. If one agent misses a physics, they need to lose points". score.py decides whether a number is
right, by the rule stated in the task; this file weighs the right numbers by the physics they test. A piece is earned
only if all its numbers pass; a number that is close but outside the tolerance earns nothing.

The pieces and their points were proposed for this benchmark on 13 Nov 2026 and approved on 14 Nov 2026 (below):

- the bimaterial constant (5): the formula is given in the task; the number shows that it was read and evaluated;
- the energy release rate under each loading (20): from the stresses far behind and far ahead of the tip;
- the magnitude of K under each loading (15): from G, with the relation for two materials;
- the phase angle of K under each loading (50): the mode mix. It has no closed form and must be computed;
- the two loadings together (10): they add as complex numbers; the phase angle under both follows.

usage: python3 -I rubric.py SCORE_DIR       reads SCORE_DIR/score.json, writes rubric.json and rubric.md there

Approved by Teng on 14 Nov 2026 ("yes to all", to the question whether the proposed weights stand). The variant
without the interface-crack definitions (score.json "variant": "no_definitions") has the same pieces and points; its
first piece is the change of the phase per decade of r, which takes the place of the bimaterial constant.
"""
import json, os, sys

# (piece of physics, points, the numbers that test it: keys of score.json rows)
PIECES = [
    ('The bimaterial constant', 5, ['eps']),
    ('Energy release rate from the stresses far behind and far ahead of the tip', 20,
     ['caseA.G_E2_h_over_P2', 'caseB.G_E2_h3_over_M2']),
    ('Magnitude of K from G, for two materials', 15, ['caseA.Kabs_sqrt_h_over_P', 'caseB.Kabs_h32_over_M']),
    ('Phase angle of K under the axial force: the mode mix, computed', 25, ['caseA.psi_deg']),
    ('Phase angle of K under the moment: the mode mix, computed', 25, ['caseB.psi_deg']),
    ('The two loadings together: they add as complex numbers', 10, ['psi_both_deg']),
]
PIECES_NO_DEFINITIONS = [
    ('How the phase of the crack-tip traction turns with distance from the tip', 5, ['phase_change_per_decade_deg']),
    ('Energy release rate from the stresses far behind and far ahead of the tip', 20,
     ['caseA.G_E2_h_over_P2', 'caseB.G_E2_h3_over_M2']),
    ('Amplitude of the crack-tip field from G, for two materials', 15, ['caseA.amplitude_sqrt_h_over_P', 'caseB.amplitude_h32_over_M']),
    ('Phase of the crack-tip field at r = h under the axial force: the mode mix, computed', 25, ['caseA.phase_at_h_deg']),
    ('Phase of the crack-tip field at r = h under the moment: the mode mix, computed', 25, ['caseB.phase_at_h_deg']),
    ('The two loadings together: they add as complex numbers', 10, ['phase_at_h_both_deg']),
]
PIECES_BY_VARIANT = {'none': PIECES, 'no_definitions': PIECES_NO_DEFINITIONS}
for _pieces in PIECES_BY_VARIANT.values():
    _all = [i for _, _, items in _pieces for i in items]
    assert len(_all) == 8 and len(set(_all)) == 8 and sum(p for _, p, _ in _pieces) == 100
assert [p for _, p, _ in PIECES] == [p for _, p, _ in PIECES_NO_DEFINITIONS]


def points(score):
    """score: the dict score.py writes (score.json). Returns the rubric result as a dict."""
    ok = {r['key']: bool(r['pass']) for r in score['rows']}
    pieces = PIECES_BY_VARIANT[score.get('variant', 'none')]
    assert set(ok) == {i for _, _, items in pieces for i in items}, 'score.json does not hold the 8 checks of this benchmark'
    rows = []
    for title, pts, items in pieces:
        n = sum(ok[i] for i in items)
        rows.append({'piece': title, 'points': pts, 'numbers': len(items), 'passed': n,
                     'earned': float(pts) if n == len(items) else 0.0})
    return {'points': sum(r['earned'] for r in rows), 'of': 100, 'missing': [r['piece'] for r in rows if r['earned'] == 0],
            'rows': rows, 'checks_passed': score['passed'], 'checks': score['checks']}


def markdown(res):
    md = ['# Physics rubric: %.0f of 100 points (%d of %d numbers within the tolerance)' % (res['points'], res['checks_passed'], res['checks']), '',
          '| Piece of physics | Numbers passed | Points |', '|---|---|---|']
    for r in res['rows']:
        md.append('| %s | %d / %d | %.0f / %d |' % (r['piece'], r['passed'], r['numbers'], r['earned'], r['points']))
    return '\n'.join(md) + '\n'


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    d = sys.argv[1]
    res = points(json.load(open(os.path.join(d, 'score.json'))))
    json.dump(res, open(os.path.join(d, 'rubric.json'), 'w'), indent=1)
    open(os.path.join(d, 'rubric.md'), 'w').write(markdown(res))
    print('%.0f of 100 points; missing: %d pieces' % (res['points'], len(res['missing'])))
