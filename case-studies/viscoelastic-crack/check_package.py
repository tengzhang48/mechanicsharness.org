"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. Runs the scorer, key/score_submission.py, on the reference values themselves (key/reference_results.json): all 251
   gates announced in task/CONTRACT.md and all 12 checks added after the audit of 25 October 2026 must pass.
4. Runs the scorer on each of the three released runs (runs/<run>/workspace/results.json) and compares its report with
   the saved one of 26 October 2026 (regrade_2026-10-26/<run>.json): every check, value and pass or fail the same. The
   explanations ("reason") are listed apart: one saved explanation predates its rewording on 26 October.
5. Runs the scorer on deliberately wrong answers made from the reference values, and confirms that exactly the expected
   gates and added checks fail. Malformed inputs must be reported, not crash the scorer.
Prints one line per check and exits with status 0 only if all pass.
"""
import copy, hashlib, json, os, shutil, subprocess, sys, tempfile

RUNS = ['P01_deepseek', 'P03_gpt61_sol', 'P04_deepseek_direct_1017']
# the saved explanations that differ from today's scorer, word for word: P04 at V = 10, which it declared inadmissible
# (FIXES.md, "Corrections after GPT's review of the push": the saved submissions score as before)
OLD_REASONS = {('P04_deepseek_direct_1017', 'V=10.0 own a, R relations'):
               ('missing finite epsilon, a or R', 'missing, non-finite or non-positive epsilon, a or R')}


def score(root, submission_path, tmp):
    """Run the scorer as a user would; return (exit status, report)."""
    report = os.path.join(tmp, 'report.json')
    if os.path.exists(report):
        os.remove(report)
    r = subprocess.run([sys.executable, '-I', os.path.join(root, 'key', 'score_submission.py'),
                        os.path.join(root, 'key', 'reference_results.json'), submission_path, report],
                       capture_output=True, text=True, timeout=120)
    if 'Traceback' in r.stderr or not os.path.exists(report):
        return r.returncode, {'scorer_crashed': (r.stderr.strip().splitlines() or ['no message'])[-1]}
    return r.returncode, json.load(open(report))


def without_reasons(report):
    out = copy.deepcopy(report)
    reasons = {}
    for key in ('checks', 'audit_checks_2026_10_26'):
        for c in out.get(key, []):
            if 'reason' in c:
                reasons[c['name']] = c.pop('reason')
    return out, reasons


def failed(report, key):
    return {c['name'] for c in report.get(key, []) if not c['passed']}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ok = True

    def report_line(passed, what):
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
        report_line(not bad and len(listed) > 0, 'SHA256SUMS: %d files present and unchanged%s' % (
            len(listed), '' if not bad else '; differ: %s' % bad[:5]))
        on_disk = {os.path.relpath(os.path.join(d, f), root).replace(os.sep, '/') for d, _, fs in os.walk(root) for f in fs}
        extra = sorted(on_disk - {p for _, p in listed} - {'SHA256SUMS'})
        report_line(not extra, 'no file outside SHA256SUMS%s' % ('' if not extra else ': %s' % extra[:5]))

        ref_path = os.path.join(root, 'key', 'reference_results.json')
        ref = json.load(open(ref_path))
        status, rep = score(root, ref_path, tmp)
        report_line(status == 0 and len(rep.get('checks', [])) == 251 and not failed(rep, 'checks')
                    and len(rep.get('audit_checks_2026_10_26', [])) == 12 and not failed(rep, 'audit_checks_2026_10_26'),
                    'the reference values pass all 251 announced gates and all 12 added checks')

        for run in RUNS:
            status, got = score(root, os.path.join(root, 'runs', run, 'workspace', 'results.json'), tmp)
            want = json.load(open(os.path.join(root, 'regrade_2026-10-26', run + '.json')))
            want = {k: v for k, v in want.items() if k not in ('run', 'announced_gates_same_as_saved')}
            g, g_reasons = without_reasons(got)
            w, w_reasons = without_reasons(want)
            reason_diff = {(run, n): (w_reasons.get(n), g_reasons.get(n)) for n in set(g_reasons) | set(w_reasons)
                           if g_reasons.get(n) != w_reasons.get(n)}
            expected_diff = {k: v for k, v in OLD_REASONS.items() if k[0] == run}
            report_line(g == w and reason_diff == expected_diff,
                        'saved answer %-25s report equals the saved one of 26 Oct (%d of %d gates and %d of %d added checks'
                        ' fail)%s' % (run, len(failed(got, 'checks')), len(got.get('checks', [])),
                                      len(failed(got, 'audit_checks_2026_10_26')), len(got.get('audit_checks_2026_10_26', [])),
                                      '; one explanation reworded since' if expected_diff else ''))

        def case(sub, v):
            return next(c for c in sub['cases'] if c['V'] == v)

        def set_value(v, key, f):
            def change(sub):
                c = case(sub, v)
                c[key] = f(c[key])
            return change

        def scale_eps_all(sub):
            for c in sub['cases']:
                c['epsilon'] *= 1.004

        def shift_sample(sub):
            case(sub, 0.1)['fields'][0]['T1'] += 0.006

        f0 = case(ref, 0.1)['fields'][0]
        sample = 'V=0.1 X=%s y/h=%s T1' % (f0['X'], f0['y_over_h'])
        # (what, change, gates that must fail, added checks that must fail)
        cases = [
            ('no bulk dissipation at V = 0.001, which the announced band passes (the audit of 25 Oct)',
             set_value(0.001, 'bulk_dissipation_ratio', lambda x: 0.0), set(),
             {'V=0.001 bulk dissipation, 5% band', 'V=0.001 energy balance, relative'}),
            ('epsilon 0.6 % high at V = 0.1, just outside the 0.5 % gate',
             set_value(0.1, 'epsilon', lambda x: 1.006 * x), {'V=0.1 epsilon'}, set()),
            ('epsilon 0.4 % high at every speed, inside the gate', scale_eps_all, set(), set()),
            ('terminal speed 1/lambda = 10, as the direct DeepSeek run inferred',
             lambda sub: sub.update(critical_speed=10.0), {'critical speed'}, set()),
            ('a 2.45 % high at V = 0.001, the size of the Azure run\'s inconsistency in its own numbers',
             set_value(0.001, 'a', lambda x: 1.0245 * x), {'V=0.001 a'},
             {'V=0.001 own a from own epsilon (Eq. 18)', 'V=0.001 own R from own a (Eq. 22)'}),
            ('negative epsilon at V = 0.1 (made the scorer crash before the fix of 26 Oct)',
             set_value(0.1, 'epsilon', lambda x: -0.5), {'V=0.1 epsilon'}, {'V=0.1 own a, R relations'}),
            ('one stress sample off by 0.006, just outside its gate', shift_sample, {sample}, set()),
        ]
        for n, (what, change, fail_gates, fail_audit) in enumerate(cases):
            sub = copy.deepcopy(ref)
            change(sub)
            path = os.path.join(tmp, 'wrong_%02d.json' % n)
            json.dump(sub, open(path, 'w'))
            status, rep = score(root, path, tmp)
            passed = (failed(rep, 'checks') == fail_gates and failed(rep, 'audit_checks_2026_10_26') == fail_audit
                      and len(rep.get('checks', [])) == 251 and status == (1 if fail_gates else 0))
            report_line(passed, ('wrong answer caught: %s (%d gates, %d added checks fail)' % (what, len(fail_gates),
                                                                                           len(fail_audit)))
                        if fail_gates or fail_audit else 'inside the gates, everything passes: %s' % what)

        path = os.path.join(tmp, 'list.json')
        json.dump([0.5447, 0.4488], open(path, 'w'))
        status, rep = score(root, path, tmp)
        report_line(status == 1 and failed(rep, 'checks') == {'top-level object'},
                    'a results file that is a list, not an object, is reported')
        path = os.path.join(tmp, 'nan.json')
        text = json.dumps(ref)
        eps = '"epsilon": %r' % case(ref, 0.1)['epsilon']
        open(path, 'w').write(text.replace(eps, '"epsilon": NaN', 1))
        status, rep = score(root, path, tmp)
        report_line(eps in text and status == 1 and 'NaN' in rep.get('format_error', ''),
                    'a results file with NaN is reported as a format error')
    print('OK: all checks passed' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
