"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. Recomputes the reference of the version without the definitions, key/reference_no_definitions.py, from
   key/reference.json and compares it with key/reference_no_definitions.json.
4. Scores each of the eight saved answers again with scorer/score.py and scorer/rubric.py (runs/<run>/submission),
   in the version the run received (runs/<run>/task_variant.txt), and compares the checks that pass and the points
   with the saved ones (runs/<run>/score).
5. Runs the scorer's controls for both versions, scorer/check_controls.py and scorer/check_controls_no_definitions.py.
Prints one line per check and exits with status 0 only if all pass.
"""
import hashlib, json, os, shutil, subprocess, sys, tempfile

RUNS = ['codex_1', 'kimi_1', 'deepseek_1', 'glm_1', 'codex_nodef1', 'kimi_nodef1', 'deepseek_nodef1', 'glm_nodef1']


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ok = True

    def report(passed, what):
        nonlocal ok
        ok &= passed
        print('%s  %s' % ('PASS' if passed else 'FAIL', what))

    def py(root, *args):
        return subprocess.run([sys.executable, '-I'] + list(args), cwd=root, capture_output=True, text=True, timeout=300)

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

        # 3. the reference of the version without the definitions
        saved = json.load(open(os.path.join(root, 'key', 'reference_no_definitions.json')))
        os.rename(os.path.join(root, 'key', 'reference_no_definitions.json'), os.path.join(tmp, 'saved.json'))
        r = py(os.path.join(root, 'key'), 'reference_no_definitions.py')
        p = os.path.join(root, 'key', 'reference_no_definitions.json')
        again = json.load(open(p)) if os.path.isfile(p) else None
        report(r.returncode == 0 and again == saved, 'key/reference_no_definitions.py rebuilds key/reference_no_definitions.json')
        shutil.copyfile(os.path.join(tmp, 'saved.json'), p)

        # 4. the saved answers, scored again
        for run in RUNS:
            out = os.path.join(tmp, 'score_' + run)
            vf = os.path.join(root, 'runs', run, 'task_variant.txt')
            variant = open(vf).read().strip() if os.path.isfile(vf) else 'none'
            r1 = py(root, os.path.join('scorer', 'score.py'), os.path.join('runs', run, 'submission'), out, variant)
            r2 = py(root, os.path.join('scorer', 'rubric.py'), out)
            if r1.returncode or r2.returncode:
                report(False, '%s: the scorer stopped: %s' % (run, (r1.stderr + r2.stderr).strip().splitlines()[-1:]))
                continue
            got, want = json.load(open(os.path.join(out, 'score.json'))), json.load(open(os.path.join(root, 'runs', run, 'score', 'score.json')))
            pg, pw = json.load(open(os.path.join(out, 'rubric.json'))), json.load(open(os.path.join(root, 'runs', run, 'score', 'rubric.json')))
            same_rows = [x['pass'] for x in got['rows']] == [x['pass'] for x in want['rows']]
            report(got['passed'] == want['passed'] and same_rows and pg['points'] == pw['points'],
                   'saved answer %-12s scores %d of %d and %.0f points, as saved' % (run, got['passed'], got['checks'], pg['points']))

        # 5. the scorer's controls, both versions
        for name in ('check_controls.py', 'check_controls_no_definitions.py'):
            r = py(os.path.join(root, 'scorer'), name)
            last = (r.stdout.strip().splitlines() or [''])[-1]
            report(r.returncode == 0 and last == 'ALL CONTROLS BEHAVE', 'scorer/%s: %s' % (name, last))

    print('OK: all checks passed' if ok else 'NOT OK: see the lines marked FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
