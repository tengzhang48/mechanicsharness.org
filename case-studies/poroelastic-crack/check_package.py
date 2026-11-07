"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. Runs the checker, key/checks.py, on each of the eighteen saved answers (runs/<run>/workspace/results.json) and compares
   its output with the saved output of 26 October 2026 (regrade_2026-10-26/<run>.json, field "checks").
4. Runs the checker on deliberately wrong answers, made by changing one value of a correct answer (MP-37), and confirms
   that exactly the expected flags change. Two changes inside the documented bands must leave every flag unchanged.
   Malformed inputs (no results file, a file that is not a JSON object) must be reported as unreadable.
Prints one line per check and exits with status 0 only if all pass.
"""
import hashlib, json, os, shutil, subprocess, sys, tempfile

RUNS = ['MP-06', 'MP-10', 'MP-13', 'MP-21', 'MP-22', 'MP-24', 'MP-27', 'MP-30', 'MP-33', 'MP-37', 'MP-38', 'MP-39',
        'MP-41', 'MP-46', 'MP-50', 'kimi-k31-high-2026-10-24', 'kimi-k31-max-2026-10-24', 'opus-5.5-subagent-2026-10-24']
K1_0, A = 0.807593, 2.440220

# (description, changes to MP-37's results.json, the flags that must change and their new values)
WRONG = [
    ('K1_0 = 1, the crack-scale value', {'K1_0': 1.0}, {'K1_0_ok': False, 'K1_0_matches': 'outer_limit_1'}),
    ('K1_0 = 0.75, the paper\'s simplified solution', {'K1_0': 0.75},
     {'K1_0_ok': False, 'K1_0_matches': 'paper_simplified_0.75'}),
    ('K1_0 3% low, outside its 2% band', {'K1_0': K1_0 * 0.97}, {'K1_0_ok': False, 'K1_0_matches': 'other'}),
    ('A = 1.354055, the face term alone', {'V1_short_coeff': 1.354055},
     {'V1_short_coeff_ok': False, 'V1_short_matches': 'faces_only_paper'}),
    ('A 6% high, outside its 5% band', {'V1_short_coeff': A * 1.06},
     {'V1_short_coeff_ok': False, 'V1_short_matches': 'other'}),
    ('K2_inf = -0.89, the effective face traction -alpha', {'K2_inf': -0.89}, {'K2_inf_ok': False}),
    ('K2_0 = 0.02, outside the absolute band 0.01', {'K2_0': 0.02}, {'K2_0_ok': False}),
    ('V2_inf with the wrong sign', {'V2_inf': 0.8388052385084748}, {'V2_inf_ok': False}),
    ('V1_0 1.5% high, outside the 1% band', {'V1_0': 1.8849555921538759 * 1.015},
     {'V1_0_ok': False, 'own_volume_ratio_is_drained_over_undrained': False}),
    ('p_faces_0 halved (the conflicting reading of the statement)', {'p_faces_0': 0.37453183520599254},
     {'p_faces_0_ok': False}),
    ('K1_inf = 0.75, not the Griffith value', {'K1_inf': 0.75}, {'K1_inf_ok': False, 'own_long_time_K_is_griffith': False}),
    ('a value given as a string', {'K2_inf': '-0.33375'}, {'K2_inf_ok': None}),
    ('K1_inf given as true, not a number', {'K1_inf': True}, {'K1_inf_ok': None, 'own_long_time_K_is_griffith': None}),
    ('K2_0 given as false, not a number', {'K2_0': False}, {'K2_0_ok': None}),
]
INSIDE = [
    ('K1_0 1.5% high, inside its 2% band', {'K1_0': K1_0 * 1.015}),
    ('A 4% low, inside its 5% band', {'V1_short_coeff': A * 0.96}),
]


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
        for run in RUNS:
            got = run_checker(root, os.path.join(root, 'runs', run, 'workspace'))
            want = json.load(open(os.path.join(root, 'regrade_2026-10-26', run + '.json')))['checks']
            n_ok = sum(1 for k, v in got.items() if k.endswith('_ok') and not k.startswith('own_') and v)
            report(got == want, 'saved answer %-29s checker output equals the saved output of 26 Oct (%d of 9 values pass)'
                   % (run, n_ok))
        base_dir = os.path.join(root, 'runs', 'MP-37', 'workspace')
        base = run_checker(root, base_dir)
        report(all(base[k] for k in base if k.endswith('_ok')), 'MP-37, the unchanged correct answer: every flag passes')
        good = json.load(open(os.path.join(base_dir, 'results.json')))
        for i, (what, change, expect) in enumerate(WRONG + [(w, c, {}) for w, c in INSIDE]):
            d = os.path.join(tmp, 'mutant_%02d' % i)
            os.makedirs(d)
            json.dump(dict(good, **change), open(os.path.join(d, 'results.json'), 'w'))
            got = run_checker(root, d)
            changed = {k: got[k] for k in got if got[k] != base[k]}
            report(changed == expect, ('wrong answer caught: %s -> %s' % (what, expect)) if expect else
                   ('inside the band, nothing changes: %s' % what))
        empty = os.path.join(tmp, 'no_results')
        os.makedirs(empty)
        got = run_checker(root, empty)
        report(list(got) == ['results_json'] and got['results_json'].startswith('unreadable'),
               'a folder without results.json is reported as unreadable')
        listfile = os.path.join(tmp, 'not_an_object')
        os.makedirs(listfile)
        json.dump([0.807593, 1.0], open(os.path.join(listfile, 'results.json'), 'w'))
        got = run_checker(root, listfile)
        report(list(got) == ['results_json'] and got['results_json'].startswith('unreadable'),
               'a results.json that is not a JSON object is reported as unreadable')
    print('OK: all checks passed' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
