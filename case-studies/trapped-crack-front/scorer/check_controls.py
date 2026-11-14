"""Controls for the scorer of the trapped-front benchmark. Exit 1 if a control does not behave as required.

1. Right answers pass every check: the key; the key rounded to 3 decimals; and the answer of this file's own solver
   (general_solve, written apart from key/reference.py) with the right coefficients.
2. Wrong theories lose checks in the groups listed for each. Every wrong theory is solved consistently by
   general_solve: its coefficients enter the stress intensity factors, the energy release rate, the equilibrium
   equation and the front together, as a solver making that mistake would have them.
3. The tolerance edge: a number moved by 0.9 of its tolerance passes, by 1.1 fails.
4. Broken files never crash the scorer and get exactly the expected number of passed checks.

general_solve. Notation (K1^0 = 1, r2 = K2^0/K1^0, r3 = K3^0/K1^0, l = L[a]/(2 pi), a' = da/dz):
  dK1 = p1 l,   dK2 = p2 r2 l - b2 r3 a',   dK3 = p3 r3 l + b3 r2 a',   G ~ K1^2 + K2^2 + m3 K3^2
  dG/G0 = (2/D) (H l + c a'),   D = 1 + r2^2 + m3 r3^2,   H = p1 + p2 r2^2 + m3 p3 r3^2,   c = r2 r3 (m3 b3 - b2).
For the right coefficients c = 0 (the slope terms cancel). Otherwise the equilibrium front also has an odd part:
with g = 2 sum q^n cos(n th), S = 2 sum q^n sin(n th), th = 2 pi s, and dG/G0 = eps g, harmonic by harmonic,
  f = (A/2pi) ln(1 - 2q cos th + q^2) + (B/pi) atan2(q sin th, 1 - q cos th),  A = D H/(H^2 + 4c^2),  B = 2 c D/(H^2 + 4c^2),
  f' = A S + B g,   l/eps = (A g - B S)/2,   and dG/(eps G0) = (2/D)(H l + c f') = g identically.

Wrong theories (named by what they change):
  meade_keer_eq16         Meade and Keer (1984b) as quoted in Gao and Rice (1986), Eq. 16: no nonlocal terms at all,
                          slope couplings 1 and 1
  rotation_scalar_kernel  shear tractions rotated into the tilted front frame (slope couplings 1 and 1) and the mode 1
                          nonlocal weight used for all three modes
  no_slope                slope coupling between modes 2 and 3 left out
  doubled_pv              nonlocal term 1/pi instead of 1/(2 pi) (the doubled gradient seen in the October pilot)
  irwin_mode3             mode 3 counted with weight 1 in G instead of 1/(1 - nu)
  slope_sign              both slope couplings with the opposite sign
  pv_swap                 the mode 2 and mode 3 nonlocal weights exchanged
Run: python3 -I check_controls.py   (standard library)
"""
import json, math, os, shutil, sys, tempfile

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE]
import score  # noqa: E402
import rubric  # noqa: E402  (points by piece of physics, on top of the checks)

REF = json.load(open(os.path.join(HERE, '..', 'key', 'reference.json')))
NU, R2, R3, Q = REF['nu'], REF['K20_over_K10'], REF['K30_over_K10'], REF['q']
RIGHT = {'p1': 1.0, 'p2': (2 - 3 * NU) / (2 - NU), 'p3': (2 + NU) / (2 - NU),
         'b2': 2 / (2 - NU), 'b3': 2 * (1 - NU) / (2 - NU), 'm3': 1 / (1 - NU)}
GA, GB1, GF, GK1, GK2, GK3, GM = [g[0] for g in score.GROUPS]     # group titles of the scorer

# name: (changed coefficients, groups in which checks must be lost, largest number of passed checks allowed)
WRONG = {
    'meade_keer_eq16': ({'p1': 0.0, 'p2': 0.0, 'p3': 0.0, 'b2': 1.0, 'b3': 1.0}, [GA, GB1, GF, GK1, GK2, GK3, GM], 20),
    'rotation_scalar_kernel': ({'p2': 1.0, 'p3': 1.0, 'b2': 1.0, 'b3': 1.0}, [GA, GB1, GF, GK1, GK2, GK3, GM], 20),
    'no_slope': ({'b2': 0.0, 'b3': 0.0}, [GA, GK2, GK3, GM], 46),
    'doubled_pv': ({'p1': 2.0, 'p2': 2 * RIGHT['p2'], 'p3': 2 * RIGHT['p3']}, [GA, GB1, GF, GK2, GK3, GM], 30),
    'irwin_mode3': ({'m3': 1.0}, [GA, GB1, GF, GK1, GK2, GK3, GM], 20),
    'slope_sign': ({'b2': -RIGHT['b2'], 'b3': -RIGHT['b3']}, [GA, GK2, GK3, GM], 46),
    'pv_swap': ({'p2': RIGHT['p3'], 'p3': RIGHT['p2']}, [GA, GB1, GF, GK1, GK2, GK3, GM], 20)}


def extremum(fun, sign=1.0, n=20000):
    """Largest value of sign * fun on -1/2 <= s < 1/2: fine scan, then golden-section refinement. Returns (s, fun(s))."""
    j = max(range(n), key=lambda i: sign * fun(i / n - 0.5))
    lo, hi = (j - 1) / n - 0.5, (j + 1) / n - 0.5
    for _ in range(200):
        a1, a2 = lo + (hi - lo) * 0.381966, lo + (hi - lo) * 0.618034
        if sign * fun(a1) < sign * fun(a2):
            lo = a1
        else:
            hi = a2
    s = 0.5 * (lo + hi)
    return s, fun(s)


def general_solve(p1, p2, p3, b2, b3, m3, r2=R2, r3=R3, q=Q):
    D = 1 + r2 ** 2 + m3 * r3 ** 2
    H = p1 + p2 * r2 ** 2 + m3 * p3 * r3 ** 2
    c = r2 * r3 * (m3 * b3 - b2)
    den = H * H + 4 * c * c
    A, B = D * H / den, 2 * c * D / den

    def dn(s):
        return 1 - 2 * q * math.cos(2 * math.pi * s) + q * q

    def g(s):
        return (1 - q * q) / dn(s) - 1

    def S(s):
        return 2 * q * math.sin(2 * math.pi * s) / dn(s)

    def f(s):
        th = 2 * math.pi * s
        return A / (2 * math.pi) * math.log(dn(s)) + B / math.pi * math.atan2(q * math.sin(th), 1 - q * math.cos(th))

    def fp(s):
        return A * S(s) + B * g(s)

    def ell(s):
        return 0.5 * (A * g(s) - B * S(s))

    def k1(s):
        return p1 * ell(s)

    def k2(s):
        return p2 * r2 * ell(s) - b2 * r3 * fp(s)

    def k3(s):
        return p3 * r3 * ell(s) + b3 * r2 * fp(s)

    # the wrong theory's own equilibrium check: dG/(eps G0) = g at a few points
    for s in (0.0, 0.13, 0.37, -0.21):
        assert abs((2 / D) * (H * ell(s) + c * fp(s)) - g(s)) < 1e-12
    st = [k / 16 for k in range(16)]
    s3, k3max = extremum(k3)
    return {'partA': {'C1': -math.pi * p1, 'S1': 0.0, 'C2': -math.pi * p2 * r2, 'S2': 2 * math.pi * b2 * r3,
                      'C3': -math.pi * p3 * r3, 'S3': -2 * math.pi * b3 * r2,
                      'CG': -2 * math.pi * H / D, 'SG': -4 * math.pi * c / D},
            'partB': {'K10_scaled': 1 / math.sqrt(D), 'front_peak_to_peak': extremum(f)[1] - extremum(f, -1.0)[1],
                      'f': [f(s) for s in st], 'kappa1': [k1(s) for s in st], 'kappa2': [k2(s) for s in st],
                      'kappa3': [k3(s) for s in st], 'kappa3_max': k3max, 's_kappa3_max': s3}}


def key_submission():
    return {'partA': dict(REF['partA']), 'partB': {k: v for k, v in REF['partB'].items() if k != 's_table'}}


def rounded(sub, nd):
    r = lambda v: [round(x, nd) for x in v] if isinstance(v, list) else round(v, nd)
    return {part: {k: r(v) for k, v in sub[part].items()} for part in sub}


def copy(sub):
    return json.loads(json.dumps(sub))


bad = 0
tmp = tempfile.mkdtemp()


def run(name, obj):
    """Score one submission. obj: a dict (written as JSON), a str (written as is) or None (empty folder)."""
    d = os.path.join(tmp, name)
    os.makedirs(d)
    if obj is not None:
        with open(os.path.join(d, 'results.json'), 'w') as fh:
            fh.write(obj if isinstance(obj, str) else json.dumps(obj))
    return score.score(d, d + '_out')


def report(ok, text):
    global bad
    bad += not ok
    print(('ok   ' if ok else 'FAIL ') + text)


try:
    # 1. right answers
    good = key_submission()
    own = general_solve(**RIGHT)
    worst = max(abs(x - y) for part in good for k in good[part]
                for x, y in zip(*[v if isinstance(v, list) else [v] for v in (good[part][k], own[part][k])]))
    report(worst < 1e-9, 'general_solve with the right coefficients against the key: largest difference %.1e' % worst)
    for name, sub in (('key', good), ('key_3_decimals', rounded(good, 3)), ('general_solve_right', own)):
        r = run(name, sub)
        report(r['passed'] == r['checks'], '%s: %d/%d' % (name, r['passed'], r['checks']))
    r = run('key_2_decimals', rounded(good, 2))
    print('note key_2_decimals: %d/%d (shows how loose the stated tolerance is; not a requirement)' % (r['passed'], r['checks']))

    # 2. wrong theories, each solved consistently
    for name, (change, must_fail, max_passed) in WRONG.items():
        r = run(name, general_solve(**dict(RIGHT, **change)))
        failing = [g['group'] for g in r['groups'] if g['passed'] < g['checks']]
        ok = r['passed'] <= max_passed and all(t in failing for t in must_fail)
        report(ok, '%s: %d/%d (at most %d allowed); ' % (name, r['passed'], r['checks'], max_passed) +
               '; '.join('%s %d/%d' % (g['group'].split(': ')[1], g['passed'], g['checks']) for g in r['groups']))

    # 3. tolerance edge: |x - ref| <= max(0.01 |ref|, 0.005)
    for part, key, idx in (('partA', 'C1', None), ('partA', 'S1', None), ('partB', 'f', 3), ('partB', 's_kappa3_max', None)):
        ref = good[part][key] if idx is None else good[part][key][idx]
        tol = max(0.01 * abs(ref), 0.005)
        for factor, want in ((0.9, 76), (1.1, 75), (-0.9, 76), (-1.1, 75)):
            s = copy(good)
            if idx is None:
                s[part][key] = ref + factor * tol
            else:
                s[part][key][idx] = ref + factor * tol
            r = run('edge_%s_%s_%+.1f' % (key, idx, factor), s)
            report(r['passed'] == want, 'edge %s%s moved by %+.1f tolerance: %d/76 (expected %d)' %
                   (key, '' if idx is None else '[%d]' % idx, factor, r['passed'], want))

    # 4. broken files: no crash, exact counts
    b = copy(good); b['partA'] = {k: False for k in b['partA']}      # false would equal the exact zeros S1 and SG
    st = copy(good); st['partB']['f'] = [str(x) for x in st['partB']['f']]
    sh = copy(good); sh['partB']['kappa2'] = sh['partB']['kappa2'][:15]
    fd = copy(good); fd['partB']['f'] = {str(k): v for k, v in enumerate(fd['partB']['f'])}
    pl = copy(good); pl['partB'] = [1, 2, 3]
    na = copy(good); del na['partA']
    ex = copy(good); ex['extra'] = 1; ex['partA']['note'] = 'x'
    huge = json.dumps(good).replace('"C1": %s' % json.dumps(good['partA']['C1']), '"C1": 1' + '0' * 400)
    nan = json.dumps(good).replace('"C1": %s' % json.dumps(good['partA']['C1']), '"C1": NaN').replace(
        '"S1": %s' % json.dumps(good['partA']['S1']), '"S1": Infinity')
    assert huge != json.dumps(good) and nan.count('NaN') == 1 and nan.count('Infinity') == 1
    broken = [('empty_folder', None, 0), ('malformed_json', '{"partA": ', 0), ('not_object', '[1, 2, 3]', 0),
              ('booleans_in_partA', b, 68), ('strings_in_f', st, 60), ('short_list_kappa2', sh, 60), ('f_is_object', fd, 60),
              ('partB_is_list', pl, 8), ('partA_missing', na, 68), ('huge_integer_C1', huge, 75),
              ('nan_and_infinity', nan, 74), ('extra_keys_ignored', ex, 76)]
    for name, obj, want in broken:
        try:
            r = run(name, obj)
            report(r['passed'] == want, '%s: %d/%d (expected %d)' % (name, r['passed'], r['checks'], want))
        except Exception as e:      # a crash of the scorer is a failed control
            report(False, '%s: the scorer crashed: %s: %s' % (name, type(e).__name__, str(e)[:80]))

    # 5. the physics rubric (rubric.py): the key earns every point; a theory that misses a piece of physics loses that
    #    piece's points and the whole outcome, however close its numbers are
    pts = lambda name, sub: rubric.points(run('rubric_' + name, sub))
    r = pts('key', good)
    report(r['points'] == 100 and r['response_complete'], 'rubric, key: %.0f of 100 points' % r['points'])
    r = pts('key_3_decimals', rounded(good, 3))
    report(r['points'] == 100, 'rubric, key rounded to 3 decimals: %.0f of 100 points' % r['points'])
    piece = {t[2][0]: t[0] for t in rubric.RESPONSE}             # first number of a piece -> its title
    # name: (pieces that must earn nothing, points expected)
    lost = {'meade_keer_eq16': (['C1', 'C2', 'C3', 'S2', 'S3', 'CG'], 10), 'rotation_scalar_kernel': (['C2', 'C3', 'S2', 'S3', 'CG'], 20),
            'no_slope': (['S2', 'S3'], 50), 'doubled_pv': (['C1', 'C2', 'C3', 'CG'], 30), 'irwin_mode3': (['K10_scaled', 'CG'], 50),
            'slope_sign': (['S2', 'S3'], 50), 'pv_swap': (['C2', 'C3', 'CG'], 40)}
    assert set(lost) == set(WRONG)
    for name, (change, _, _) in WRONG.items():
        r = pts(name, general_solve(**dict(RIGHT, **change)))
        earned = {x['piece']: x['earned'] for x in r['rows']}
        out = sum(x['earned'] for x in r['rows'] if x['part'] == 'outcome')
        ok = (all(earned[piece[k]] == 0 for k in lost[name][0]) and out == 0 and r['points'] == lost[name][1]
              and sum(v == 0 for t, v in earned.items() if t in piece.values()) == len(lost[name][0]))
        report(ok, 'rubric, %s: %.0f of 100 points (expected %d), outcome %.0f, %d pieces of the response lost' %
               (name, r['points'], lost[name][1], out, len(r['missing'])))
    # a slip in one number of a table, with the whole response right, costs a fraction of a point, not the outcome
    slip = copy(good)
    slip['partB']['kappa2'][5] += 1.0
    r = pts('one_table_value_wrong', slip)
    report(r['response_complete'] and abs(r['points'] - (100 - 20 / 50)) < 1e-9, 'rubric, one table value wrong: %.1f of 100 points' % r['points'])
finally:
    shutil.rmtree(tmp)
print('FAILED %d' % bad if bad else 'ALL CONTROLS BEHAVE')
sys.exit(1 if bad else 0)
