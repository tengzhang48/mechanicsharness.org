"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. Recomputes the sixteen key values from the closed forms in key/reference_frozen/REFERENCE.md and compares them with the
   frozen reference, key/reference_frozen/reference.json.
4. Runs the checker, key/checks.py, on each of the five released runs (runs/<run>/workspace) and compares its output with
   the saved output of 26 October 2026 (regrade_2026-10-26/<run>.json).
5. Runs the checker on deliberately wrong answers, those the audit of 25 October 2026 tried (audit_2026-10-25/audit.json),
   and confirms that exactly the expected values fail. A change inside the 0.5 % band must pass every value, and one just
   outside it must fail every value. Values that are not numbers, and malformed inputs, must be reported.
Prints one line per check and exits with status 0 only if all pass.
"""
import hashlib, json, math, os, shutil, subprocess, sys, tempfile

RUNS = ['S01', 'S02', 'S04', 'D01', 'K01-continuation']
NU = 0.3
LH = (0.5, 1.0, 1.5, 2.0)   # the sampling locations ell/H
Q = ('gamma_A', 'gamma_B', 'kI_A', 'kI_B')


def closed_forms(x):
    """x = H/ell. gamma = G Ebar/(sigma0^2 ell); kI = KI/(sigma0 sqrt(ell)) = sqrt((1 - nu^2) gamma), since G = KI^2/E."""
    gA, gB = math.tanh(x / 2), 2 * math.tanh(x / 2) - math.tanh(x)
    return {'gamma_A': gA, 'gamma_B': gB, 'kI_A': math.sqrt((1 - NU ** 2) * gA), 'kI_B': math.sqrt((1 - NU ** 2) * gB)}


def flags(qs, lhs=LH):
    return {'%s_%s_ok' % (q, lh) for q in qs for lh in lhs}


def run_checker(root, folder):
    r = subprocess.run([sys.executable, '-I', os.path.join(root, 'key', 'checks.py'), folder], capture_output=True,
                       text=True, timeout=60)
    if r.returncode:  # a crash is a failed check, reported with its last line
        return {'checker_crashed': (r.stderr.strip().splitlines() or ['no message'])[-1]}
    return json.loads(r.stdout)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ok = True

    def report(passed, what):
        nonlocal ok
        ok &= passed
        print('%s  %s' % ('PASS' if passed else 'FAIL', what))

    with tempfile.TemporaryDirectory() as tmp:
        root = os.path.join(tmp, 'package')
        shutil.copytree(here, root)
        print('working in a clean copy: %s' % root)
        listed = [line.split('  ', 1) for line in open(os.path.join(root, 'SHA256SUMS')).read().splitlines()]
        bad = [p for h, p in listed if not os.path.isfile(os.path.join(root, p)) or
               hashlib.sha256(open(os.path.join(root, p), 'rb').read()).hexdigest() != h]
        report(not bad and len(listed) > 0, 'SHA256SUMS: %d files present and unchanged%s' % (
            len(listed), '' if not bad else '; differ: %s' % bad[:5]))
        on_disk = {os.path.relpath(os.path.join(d, f), root).replace(os.sep, '/') for d, _, fs in os.walk(root) for f in fs}
        extra = sorted(on_disk - {p for _, p in listed} - {'SHA256SUMS'})
        report(not extra, 'no file outside SHA256SUMS%s' % ('' if not extra else ': %s' % extra[:5]))

        ref = json.load(open(os.path.join(root, 'key', 'reference_frozen', 'reference.json')))['rows']
        key = [dict(ell_over_H=lh, **closed_forms(1 / lh)) for lh in LH]
        dev = max(abs(k[q] - r[q]) for k, r in zip(key, ref) for q in Q)
        report([r['ell_over_H'] for r in ref] == list(LH) and dev < 1e-12,
               'the closed forms give the frozen reference values (largest difference %.1e)' % dev)

        for run in RUNS:
            got = run_checker(root, os.path.join(root, 'runs', run, 'workspace'))
            want = json.load(open(os.path.join(root, 'regrade_2026-10-26', run + '.json')))
            report(got == want, 'saved answer %-17s checker output equals the saved output of 26 Oct (%s of 16 within 0.5 %%)'
                   % (run, got.get('values_within_0.5pct')))

        controls = json.load(open(os.path.join(root, 'key', 'reference_frozen', 'reference.json')))['meta'][
            'controls_deliberately_wrong']

        def rows(change):
            out = []
            for i, k in enumerate(key):
                r = dict(k)
                change(i, r)
                out.append(r)
            return out

        def with_gamma(name, values):  # a wrong gamma and the kI that follows from it
            def change(i, r):
                r['gamma_' + name] = values[i]
                r['kI_' + name] = math.sqrt((1 - NU ** 2) * values[i])
            return change

        def at(xfun):  # every value at the wrong argument
            return lambda i, r: r.update(closed_forms(xfun(LH[i])))

        # (description, change, values that must fail, values that must be reported missing)
        cases = [
            ('kI from G with Ebar, without 1 - nu^2',
             lambda i, r: r.update(kI_A=r['kI_A'] / math.sqrt(1 - NU ** 2), kI_B=r['kI_B'] / math.sqrt(1 - NU ** 2)),
             flags(('kI_A', 'kI_B')), []),
            ('foundation energy left out (the reference\'s control C1)',
             with_gamma('A', controls['C1_omit_foundation_energy_gamma_A_wrong']), flags(('gamma_A', 'kI_A')), []),
            ('G of case B per period, not per tip (x 2)',
             lambda i, r: r.update(gamma_B=2 * r['gamma_B'], kI_B=math.sqrt(2) * r['kI_B']), flags(('gamma_B', 'kI_B')), []),
            ('case B spacings swapped (the reference\'s control C3)',
             with_gamma('B', controls['C3_swapped_spacings_gamma_B_wrong']), flags(('gamma_B', 'kI_B')), []),
            ('ell/H read as H/ell', at(lambda lh: lh), flags(Q, (0.5, 1.5, 2.0)), []),
            ('decay length 2 ell (H/4ell in place of H/2ell)', at(lambda lh: 1 / (2 * lh)), flags(Q), []),
            ('every value 0.6 % high, just outside the band', lambda i, r: r.update({q: 1.006 * r[q] for q in Q}), flags(Q), []),
            ('gamma_A at ell/H = 0.5 given as true, not a number', lambda i, r: r.update(gamma_A=True) if i == 0 else None,
             flags(('gamma_A',), (0.5,)), ['gamma_A at ell/H = 0.5']),
            ('ell/H = 1.5 given as a string', lambda i, r: r.update(ell_over_H='1.5') if i == 2 else None,
             flags(Q, (1.5,)), ['%s at ell/H = 1.5' % q for q in Q]),
            ('every value 0.4 % high, inside the band', lambda i, r: r.update({q: 1.004 * r[q] for q in Q}), set(), []),
        ]
        d = os.path.join(tmp, 'key_rows')
        os.makedirs(d)
        json.dump({'rows': key}, open(os.path.join(d, 'results.json'), 'w'))
        got = run_checker(root, d)
        report(got.get('values_within_0.5pct') == 16 and got.get('values_within_1e-6_abs') == 16 and not got.get('missing'),
               'the closed-form values: all 16 pass, also within 1e-6')
        for n, (what, change, fail, missing) in enumerate(cases):
            d = os.path.join(tmp, 'wrong_%02d' % n)
            os.makedirs(d)
            json.dump({'rows': rows(change)}, open(os.path.join(d, 'results.json'), 'w'))
            got = run_checker(root, d)
            failed = {k for k, v in got.items() if k.endswith('_ok') and v is False}
            passed = (failed == fail and got.get('values_within_0.5pct') == 16 - len(fail) and got.get('missing') == missing
                      and len([k for k in got if k.endswith('_ok')]) == 16)
            report(passed, ('wrong answer caught: %s (%d of 16 fail)' % (what, len(fail))) if fail else
                   ('inside the band, every value passes: %s' % what))
        empty = os.path.join(tmp, 'no_results')
        os.makedirs(empty)
        got = run_checker(root, empty)
        report(list(got) == ['results_json'], 'a folder without results.json or numeric_rows.json is reported')
        numbers = os.path.join(tmp, 'numbers_only')
        os.makedirs(numbers)
        json.dump([0.7615941559557649, 0.5591607318357128], open(os.path.join(numbers, 'results.json'), 'w'))
        got = run_checker(root, numbers)
        report(got.get('values_within_0.5pct') == 0 and len(got.get('missing', [])) == 16,
               'a results.json that is a list of numbers, not rows: all 16 values reported missing')
    print('OK: all checks passed' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
