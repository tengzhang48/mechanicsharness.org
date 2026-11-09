"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. Reruns the checker's self-test (key/check_P038.py --selftest: the key, the audit's independent recomputation, the
   numbers of the eight graded answers and the audit's mutants) and compares its output with the saved one
   (regrade_2026-10-26/check_P038_selftest.json).
4. Runs the checker on the key, on values just inside the full-credit bands, on deliberately wrong answers, on a value
   given as text and on inputs with no problem block, and confirms what it reports.
Prints one line per check and exits with status 0 only if all pass.
"""
import hashlib, json, os, shutil, subprocess, sys, tempfile

KEY = {'T1': {'Sigma_MPa': 133.730286053811, 'L_mm_per_end': 4.08560913824553},
       'T2': {'M_Nm_each_strip': 0.212625787376721, 'L_mm': 10.8001023708973, 'delta_mid_um': 38.2188979260933}}


def run(args, cwd):
    return subprocess.run([sys.executable, '-I'] + args, cwd=cwd, capture_output=True, text=True, timeout=120)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ok = True

    def report(passed, what):
        nonlocal ok
        ok &= bool(passed)
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

        out = os.path.join(tmp, 'selftest.json')
        r = run(['key/check_P038.py', '--selftest', '--out', out], root)
        saved = json.load(open(os.path.join(root, 'regrade_2026-10-26', 'check_P038_selftest.json')))
        got = json.load(open(out)) if r.returncode == 0 and os.path.exists(out) else None
        report(got == saved, 'checker self-test (%d cases: key, independent recomputation, eight graded answers, %d mutants)'
               ' equals the saved one' % (len(saved), sum('mutant' in c['name'] for c in saved)))
        if got:
            pts = {c['name']: c['points'] for c in got}
            report(all(v == {'T1': 50} for k, v in pts.items() if k.startswith('run:P038-T1') and 'R03' not in k)
                   and pts['run:P038-T1-without-dsaz (R03)'] == {'T1': 0}
                   and all(v == {'T2': 40} for k, v in pts.items() if k.startswith('run:P038-T2')),
                   'saved answers: T1 R08, R24, R30 50 of 50 and R03 (wrong opening) 0 of 50; all four T2 answers 40 of 40')
            report(not any(all(all(f.values()) for f in c['pass_flags'].values()) for c in got if 'mutant' in c['name']),
                   'every one of the audit\'s mutants fails at least one full-credit band')

        def check(values):
            f = os.path.join(tmp, 'r_%d.json' % len(os.listdir(tmp)))
            json.dump(values, open(f, 'w'))
            r = run(['key/check_P038.py', f], root)
            return r.returncode, (json.loads(r.stdout) if r.returncode == 0 else r.stdout.strip())

        code, got = check(KEY)
        report(code == 0 and got['T1']['numeric_points_of_50'] == 50 and got['T2']['numeric_points_of_40'] == 40,
               'the key: 50 of 50 (T1) and 40 of 40 (T2)')
        inside = {'T1': {k: v * 1.0025 for k, v in KEY['T1'].items()}, 'T2': {k: v * 1.004 for k, v in KEY['T2'].items()}}
        code, got = check(inside)
        report(code == 0 and got['T1']['numeric_points_of_50'] == 50 and got['T2']['numeric_points_of_40'] == 40,
               'values just inside the full-credit bands (T1 0.25 % high, T2 0.4 % high) get full points')
        for what, values in [('T1: R03\'s numbers (158.5 MPa, 34.3 mm)', {'T1': {'Sigma_MPa': 158.5, 'L_mm_per_end': 34.3}}),
                             ('T1: Sigma 1 % high', {'T1': dict(KEY['T1'], Sigma_MPa=KEY['T1']['Sigma_MPa'] * 1.01)}),
                             ('T2: moment summed over both strips', {'T2': dict(KEY['T2'], M_Nm_each_strip=2 * KEY['T2']['M_Nm_each_strip'])})]:
            code, got = check(values)
            prob = list(values)[0]
            report(code == 0 and not all(v['pass'] for k, v in got[prob].items() if isinstance(v, dict)),
                   'wrong answer caught: %s' % what)
        code, got = check({'T1': dict(KEY['T1'], Sigma_MPa='133.73')})
        report(code == 0 and got['T1']['Sigma_MPa'].get('note') == 'not a number' and got['T1']['L_mm_per_end']['pass'],
               'a value given as text is reported as not a number, the other still checked')
        for values in ({}, [], {'T9': {'x': 1.0}}):
            code, got = check(values)
            report(code == 1 and 'nothing was checked' in got, 'an input with no problem block (%s) is reported, not passed'
                   % json.dumps(values))
    print('OK: all checks passed' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
