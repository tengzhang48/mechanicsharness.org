"""Check this case package from a clean copy. Python 3.8 or later, standard library only.

    python3 check_package.py

1. Copies this folder to a new temporary directory and works only there.
2. Checks every file listed in SHA256SUMS: present and unchanged, and no file that is not listed.
3. T1 and T2: reruns the checker's self-test (key/check_P015.py --selftest: the key, the numbers of the six graded
   answers and the audit's mutants) and compares its output with the saved one (regrade_2026-10-26/check_P015_selftest.json).
4. V1 and V2: reruns the checker's audit (matched_test/key/check_V1_V2.py --audit, with the audit's transcribed answers and
   mutants) and compares its output with the saved one (matched_test/audit_2026-11-07/check_V1_V2_audit.out).
5. Runs both checkers on the keys, on values just inside the bands, on deliberately wrong answers, on a value given as
   text and on inputs with no problem block, and confirms what each reports.
Prints one line per check and exits with status 0 only if all pass.
"""
import hashlib, json, os, shutil, subprocess, sys, tempfile


def run(args, cwd):
    return subprocess.run([sys.executable, '-I'] + args, cwd=cwd, capture_output=True, text=True, timeout=120)


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

        # T1 and T2: the checker's self-test against the saved one
        out = os.path.join(tmp, 'selftest.json')
        r = run(['key/check_P015.py', '--selftest', '--out', out], root)
        saved = json.load(open(os.path.join(root, 'regrade_2026-10-26', 'check_P015_selftest.json')))
        got = json.load(open(out)) if r.returncode == 0 and os.path.exists(out) else None
        report(got == saved, 'T1/T2 checker self-test (%d cases: key, six graded answers, %d mutants) equals the saved one'
               % (len(saved), sum(c['name'].startswith('mutant') for c in saved)))
        if got:
            graded = {c['name']: c['checker_result'] for c in got if c['name'].startswith('run ')}
            report(all(v[t]['all_full'] for k, v in graded.items() if 'R11' in k or 'R35' in k or 'R32' in k or 'R09' in k
                       for t in v) and not any(v[t]['all_full'] for k, v in graded.items()
                                               if 'R31' in k or 'R23' in k for t in v),
                   'saved answers: R11, R35, R32 and R09 at full credit; R31 (J without the shear-traction work) and R23'
                   ' (wrong notch equations) not')
            report(not any(c['checker_result'][t]['all_full'] for c in got if c['name'].startswith('mutant')
                           for t in c['checker_result']), 'every one of the audit\'s T1/T2 mutants is caught')

        def p015(values):
            f = os.path.join(tmp, 'p015_%d.json' % len(os.listdir(tmp)))
            json.dump(values, open(f, 'w'))
            r = run(['key/check_P015.py', f], root)
            return json.loads(r.stdout) if r.returncode == 0 else {'crashed': r.stderr.strip()[-200:]}

        key = {'T1': {'p': 2 / 15, 'H': 1.893569725213, 'C': 0.457305783974}, 'T2': {'p': 0.156175260286, 'R': 0.373262313312}}
        got = p015(key)
        report(got.get('T1', {}).get('all_full') and got.get('T2', {}).get('all_full'), 'T1/T2 keys: every number at full credit')
        got = p015({'T1': {'p': 2 / 15, 'H': 1.8975, 'C': 0.4588}, 'T2': {'p': 0.15645, 'R': 0.3750}})
        report(got.get('T1', {}).get('all_full') and got.get('T2', {}).get('all_full'),
               'T1/T2 values just inside the full-credit bands pass')
        mut = json.load(open(os.path.join(root, 'audit_2026-10-25', 'mutant_values.json')))
        for task, name, what in [('T1', 'plane_stress_approx', 'plane stress instead of plane strain'),
                                 ('T1', 'C_without_n_over_n1', 'C without the factor n/(n+1)'),
                                 ('T2', 'crack_exponent_n4', 'the crack exponent 1/(n+1) at the notch'),
                                 ('T2', 'elastic_wedge', 'the elastic wedge exponent at the notch')]:
            v = {k: mut[task][name][k] for k in (('p', 'H', 'C') if task == 'T1' else ('p', 'R'))}
            got = p015({task: v})
            report(got.get(task, {}).get('all_full') is False, 'wrong answer caught: %s %s' % (task, what))
        got = p015({'T1': {'p': '0.1333', 'H': 1.893569725213, 'C': 0.457305783974}})
        report(got.get('T1', {}).get('p') == 'not a number' and got['T1']['all_full'] is False and got['T1']['H'] == 'full',
               'T1/T2: a value given as text is reported as not a number, the others still checked')

        # V1 and V2: the checker's audit against the saved one
        work = os.path.join(tmp, 'v1v2')
        os.makedirs(work)
        for f in ('key/check_V1_V2.py', 'audit_2026-11-07/answers_extracted.json', 'audit_2026-11-07/mutants.json'):
            shutil.copy(os.path.join(root, 'matched_test', f), work)
        r = run(['check_V1_V2.py', '--audit'], work)
        saved = open(os.path.join(root, 'matched_test', 'audit_2026-11-07', 'check_V1_V2_audit.out')).read()
        got = open(os.path.join(work, 'check_V1_V2_audit.out')).read() if r.returncode == 0 else ''
        report(got == saved, 'V1/V2 checker audit (key, 24 graded answers, 3 Kimi answers, mutants) equals the saved one')
        agree = [l for l in got.splitlines() if ' checker ' in l and ' grader ' in l and ' M' in l.split('P015')[0]]
        report(len(agree) == 24 and all(l.split('checker ')[1][:3] == l.split('grader ')[1][:3] for l in agree),
               'V1/V2: the checker agrees with all 24 graders\' tolerance verdicts')

        def v1v2(values):
            f = os.path.join(tmp, 'v_%d.json' % len(os.listdir(tmp)))
            json.dump(values, open(f, 'w'))
            r = run(['matched_test/key/check_V1_V2.py', f], root)
            fails = [l.split()[1] for l in r.stdout.splitlines() if l.endswith('FAIL')]
            return r.returncode, fails

        kv = json.load(open(os.path.join(root, 'matched_test', 'audit_2026-11-07', 'key_from_recompute.json')))
        code, fails = v1v2(kv)
        report(code == 0 and not fails, 'V1/V2 keys as recomputed by the audit: all fourteen numbers pass')
        inside = {'V1': {k: v * 1.015 for k, v in kv['V1'].items()}, 'V2': dict(kv['V2'], T=kv['V2']['T'] * 1.008)}
        code, fails = v1v2(inside)
        report(code == 0, 'V1/V2 values just inside their bands (V1 1.5 % high, V2 T 0.8 % high) pass')
        mv = json.load(open(os.path.join(root, 'matched_test', 'audit_2026-11-07', 'mutants.json')))
        for name, what in [('V1_law_applied_to_tensor_shear', 'the law applied to tensor shear'),
                           ('V1_J_term_dropped_traction', 'the traction term dropped from J'),
                           ('V2_incompressible_nu_t_half_paper_eq48', 'the incompressible approximation'),
                           ('V2_plane_strain_J_K_relation', 'the plane-strain J-K relation')]:
            prob = name[:2]
            values = mv[name]['values'] if isinstance(mv[name], dict) and 'values' in mv[name] else mv[name]
            code, fails = v1v2({prob: values})
            report(code == 1 and bool(fails), 'wrong answer caught: %s %s (%s fail)' % (prob, what, ', '.join(fails)))
        code, fails = v1v2({'V1': dict(kv['V1'], tau_rz='38.59')})
        report(code == 1 and fails == ['tau_rz'], 'V1/V2: a value given as text fails, the others still pass')
        # inputs with no problem block: both checkers must say so and exit with status 1
        for values in ({}, [], {'V3': {'x': 1.0}}):
            f = os.path.join(tmp, 'empty_%d.json' % len(os.listdir(tmp)))
            json.dump(values, open(f, 'w'))
            a = run(['key/check_P015.py', f], root)
            b = run(['matched_test/key/check_V1_V2.py', f], root)
            report(a.returncode == 1 and b.returncode == 1 and 'nothing was checked' in a.stdout + b.stdout,
                   'an input with no problem block (%s) is reported, not passed, by both checkers' % json.dumps(values))
    print('OK: all checks passed' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
