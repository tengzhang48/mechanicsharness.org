"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. Runs the checker, key/checks.py, on each of the eleven runs it was run on (runs/<run>/workspace) and compares its flags
   with the saved ones of 26 October 2026 (regrade_2026-10-26/<run>.json). For the runs that saved no results.json the
   saved output names a temporary folder; there the check is that both report the file as unreadable.
4. Runs the checker on the key's own values (all flags pass, lambda1 "exact") and on deliberately wrong answers,
   including those of FIXES.md (the paper's Appendix B estimates lambda1 = 13.36 and shape 0.811, and K = 1), and
   confirms that exactly the expected flags fail.
Prints one line per check and exits with status 0 only if all pass.
"""
import hashlib, json, os, shutil, subprocess, sys, tempfile

RUNS = ['MP-02', 'MP-03', 'MP-05', 'MP-09', 'MP-15', 'MP-18', 'MP-20', 'MP-32', 'MP-40', 'MP-44', 'MP-47']
KEY = {'lambda1': 11.742662, 'lambda2': 264.6346, 'mode1_shape_half': 0.80196, 'K_inf_norm': 1.121522,
       'relaxation_time_t0': 0.0852}   # key/KEY.md: lambda1, lambda2, f1(0.5), the traction-free limit, the e-folding time


def check(root, folder):
    r = subprocess.run([sys.executable, '-I', os.path.join(root, 'key', 'checks.py'), folder], capture_output=True,
                       text=True, timeout=60)
    if r.returncode or 'Traceback' in r.stderr:
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

        for run in RUNS:
            got = check(root, os.path.join(root, 'runs', run, 'workspace'))
            want = json.load(open(os.path.join(root, 'regrade_2026-10-26', run + '.json')))['checks']
            if 'results_json' in want:
                same = list(got) == ['results_json'] and got['results_json'].startswith('unreadable')
                what = 'no results.json, reported as unreadable'
            else:
                same = got == want
                what = 'lambda1 %s, %d of %d flags pass' % (got.get('lambda1_class'), sum(v is True for v in got.values()),
                                                            sum(isinstance(v, bool) for v in got.values()))
            report(same, 'saved answer %-6s checker output equals the saved one of 26 Oct (%s)' % (run, what))

        def flags(values):
            d = os.path.join(tmp, 'answer_%d' % len(os.listdir(tmp)))
            os.makedirs(d)
            if values is not None:
                json.dump(values, open(os.path.join(d, 'results.json'), 'w'))
            return check(root, d)

        got = flags(KEY)
        report(got.get('lambda1_class') == 'exact' and all(v is True for k, v in got.items() if k != 'lambda1_class'),
               'the key\'s own values: lambda1 exact and every flag passes')
        # (what, change, flags that must be false, lambda1 class)
        cases = [
            ('lambda1 = 13.36, the paper\'s crack-shape estimate (Appendix B)', {'lambda1': 13.36}, {'lambda1_ok'}, 'approximate'),
            ('mode shape 0.811, the paper\'s Appendix B estimate', {'mode1_shape_half': 0.811}, {'mode1_shape_ok'}, 'exact'),
            ('K = 1, without the edge-crack factor 1.1215', {'K_inf_norm': 1.0}, {'K_inf_ok'}, 'exact'),
            ('lambda1 1.2 % high, just outside the 1 % band', {'lambda1': 1.012 * KEY['lambda1']}, {'lambda1_ok'}, 'approximate'),
            ('lambda2 below lambda1', {'lambda2': 5.0}, {'lambda2_ok', 'own_lambda2_gt_lambda1'}, 'exact'),
            ('relaxation time 0.3, outside the 0.05-0.15 band and inconsistent with lambda1', {'relaxation_time_t0': 0.3},
             {'relaxation_time_ok', 'own_time_consistent_with_lambda1'}, 'exact'),
            ('lambda1 0.6 % high, inside the 1 % band', {'lambda1': 1.006 * KEY['lambda1']}, set(), 'exact'),
        ]
        for what, change, fail, cls in cases:
            got = flags(dict(KEY, **change))
            failed = {k for k, v in got.items() if v is False}
            report(failed == fail and got.get('lambda1_class') == cls and all(v is not None for v in got.values()),
                   ('wrong answer caught: %s (%s)' % (what, ', '.join(sorted(fail)))) if fail else
                   'inside the bands, every flag passes: %s' % what)
        got = flags(dict(KEY, lambda2='264.6'))
        report(got.get('lambda2_ok') is None and got.get('own_lambda2_gt_lambda1') is None and got.get('lambda1_ok') is True,
               'a value given as text is not counted (its flags are null), the others still are')
        got = flags(None)
        report(list(got) == ['results_json'], 'a folder without results.json is reported as unreadable')
    print('OK: all checks passed' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
