"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. Recomputes the reference, key/reference.py, and compares it with key/reference.json.
4. Scores each of the seven saved answers again with scorer/score.py and scorer/rubric.py (runs/<run>/submission) and
   compares the checks that pass and the points with the saved ones (runs/<run>/score). The run glm_r1 saved no
   results.json: there the check is that it scores 0 of 76 and 0 points, as saved.
5. Runs the scorer's controls, scorer/check_controls.py: the key earns 76 of 76 and 100 points, seven wrong first-order
   theories lose the checks and the pieces of physics they miss, the tolerance sits where the task states it, and broken
   files fail without stopping the scorer.
Prints one line per check and exits with status 0 only if all pass.
"""
import hashlib, json, os, shutil, subprocess, sys, tempfile

RUNS = ['codex_1', 'kimi_1', 'deepseek_1', 'glm_1', 'kimi_r1', 'deepseek_r1', 'glm_r1']


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

        # 3. the reference
        saved = json.load(open(os.path.join(root, 'key', 'reference.json')))
        os.rename(os.path.join(root, 'key', 'reference.json'), os.path.join(tmp, 'reference.saved.json'))
        r = py(os.path.join(root, 'key'), 'reference.py')
        again = json.load(open(os.path.join(root, 'key', 'reference.json'))) if os.path.isfile(os.path.join(root, 'key', 'reference.json')) else None

        def flat(x, pre=''):
            if isinstance(x, dict):
                return {k2: v2 for k, v in x.items() for k2, v2 in flat(v, pre + k + '.').items()}
            if isinstance(x, list):
                return {k2: v2 for i, v in enumerate(x) for k2, v2 in flat(v, pre + str(i) + '.').items()}
            return {pre.rstrip('.'): x}
        fs, fa = flat(saved), flat(again) if again else {}
        num = [k for k, v in fs.items() if isinstance(v, (int, float)) and not isinstance(v, bool)]
        worst = max((abs(fa[k] - fs[k]) for k in num if k in fa), default=float('inf'))
        report(r.returncode == 0 and again is not None and set(fs) == set(fa) and worst < 1e-12,
               'key/reference.py rebuilds key/reference.json (%d numbers, largest difference %.1e)' % (len(num), worst))
        shutil.copyfile(os.path.join(tmp, 'reference.saved.json'), os.path.join(root, 'key', 'reference.json'))

        # 4. the saved answers, scored again
        for run in RUNS:
            out = os.path.join(tmp, 'score_' + run)
            r1 = py(root, os.path.join('scorer', 'score.py'), os.path.join('runs', run, 'submission'), out)
            r2 = py(root, os.path.join('scorer', 'rubric.py'), out)
            if r1.returncode or r2.returncode:
                report(False, '%s: the scorer stopped: %s' % (run, (r1.stderr + r2.stderr).strip().splitlines()[-1:]))
                continue
            got, want = json.load(open(os.path.join(out, 'score.json'))), json.load(open(os.path.join(root, 'runs', run, 'score', 'score.json')))
            pg, pw = json.load(open(os.path.join(out, 'rubric.json'))), json.load(open(os.path.join(root, 'runs', run, 'score', 'rubric.json')))
            same_rows = [x['pass'] for x in got['rows']] == [x['pass'] for x in want['rows']]
            report(got['passed'] == want['passed'] and same_rows and pg['points'] == pw['points'],
                   'saved answer %-12s scores %d of %d and %.0f points, as saved' % (run, got['passed'], got['checks'], pg['points']))

        # 5. the scorer's controls
        r = py(os.path.join(root, 'scorer'), 'check_controls.py')
        last = (r.stdout.strip().splitlines() or [''])[-1]
        report(r.returncode == 0 and last == 'ALL CONTROLS BEHAVE', 'scorer/check_controls.py: %s' % last)

    print('OK: all checks passed' if ok else 'NOT OK: see the lines marked FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
