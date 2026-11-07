"""Master driver: runs all checks and assembles results.json.

Single reproduction command:  python3 run_all.py
"""
import json
import math
import subprocess
import sys
import re

NU = 0.3


def run(cmd):
    p = subprocess.run([sys.executable, cmd], capture_output=True, text=True, timeout=600)
    return p.stdout + p.stderr


def main():
    # ---- analytic sampling rows -------------------------------------------------
    rows = []
    for eoH in (0.5, 1.0, 1.5, 2.0):
        x = 1.0 / (2.0 * eoH)
        t = math.tanh(x)
        gA = t
        gB = 2.0 * t - math.tanh(2.0 * x)
        rows.append(dict(ell_over_H=eoH, H_over_2ell=x,
                         gamma_A=gA, gamma_B=gB,
                         kI_A=math.sqrt((1 - NU**2) * gA),
                         kI_B=math.sqrt((1 - NU**2) * gB)))

    checks = []

    # ---- 1. symbolic checks ------------------------------------------------------
    out = run('analytic_checks.py')
    open('analytic_output.txt', 'w').write(out)
    checks.append(dict(
        name='sympy_symbolic_verification',
        description='wake ODE residual=0; BC sigma22=0 at both strip edges; strip energy '
                    'integral W(S) = -h s0^2 ell tanh(S/2ell)/Ebar (exact + 4-point numeric); '
                    'linear/quadratic decomposition -2:+1; identity 2tanh x - tanh 2x = 2t^3/(1+t^2)',
        output='8/8 PASS' if 'ALL SYMBOLIC CHECKS PASS' in out else 'FAILURE: ' + out[-500:],
        tolerance='exact (0 residual)', passed='ALL SYMBOLIC CHECKS PASS' in out))

    # ---- 2. numeric eval (identity of the two gamma_B forms) ---------------------
    out = run('numeric_eval.py')
    open('numeric_output.txt', 'w').write(out)
    checks.append(dict(
        name='gamma_B_two_forms_consistency',
        description='2 tanh(x) - tanh(2x) vs 2 t^3/(1+t^2) at the four sampling points',
        output='max |diff| < 1e-15 (assert in numeric_eval.py) ' + ('OK' if 'wrote' in out else 'FAIL'),
        tolerance=1e-15, passed='wrote' in out))

    # ---- 3. FD J-integral check ----------------------------------------------------
    out = run('fd_check.py')
    open('fd_output.txt', 'w').write(out)
    fd_rows = re.findall(r'case (\w) H=\s*([\d.]+).*?rel_err=([+-][\d.]+)%', out)
    conv = re.search(r'ratio=([\d.]+)', out)
    fd_pass = 'FD CHECKS PASS' in out
    checks.append(dict(
        name='finite_difference_J_integral',
        description='independent 2D FD solution on symmetry half-domain; G = (2/h) DJ '
                    'from path-independent Eshelby line integrals; 5 configurations',
        output='; '.join(f'case {c} H={h}: rel_err {e}%' for c, h, e in fd_rows)
               + (f'; grid-convergence ratio {conv.group(1)} (expect ~4 for 2nd order)' if conv else ''),
        tolerance='|rel err| < 3%', passed=fd_pass))

    # ---- 4. asymptotic limits --------------------------------------------------------
    xs = [1e-1, 1e-2, 1e-3]
    smallA = [math.tanh(x) / x for x in xs]                       # -> 1
    smallB = [(2 * math.tanh(x) - math.tanh(2 * x)) / (2 * x**3) for x in xs]  # -> 1
    big = 20.0
    bigA = math.tanh(big)                                          # -> 1
    bigB = 2 * math.tanh(big) - math.tanh(2 * big)                 # -> 1
    lim_ok = (abs(smallA[-1] - 1) < 1e-5 and abs(smallB[-1] - 1) < 1e-4
              and abs(bigA - 1) < 1e-15 and abs(bigB - 1) < 1e-15)
    checks.append(dict(
        name='asymptotic_limits',
        description='gamma_A/x -> 1 and gamma_B/(2x^3) -> 1 as x=H/2ell -> 0; '
                    'gamma_A, gamma_B -> 1 as x -> infinity',
        output=(f'gamma_A/x at x={xs}: {[round(v, 8) for v in smallA]}; '
                f'gamma_B/(2x^3) at x={xs}: {[round(v, 8) for v in smallB]}; '
                f'gamma_A(x=20)={bigA:.16f}, gamma_B(x=20)={bigB:.16f}'),
        tolerance='ratios -> 1 within 1e-4 at x=1e-3; big-x within 1e-15', passed=lim_ok))

    results = dict(
        problem='Energy release for periodic cracks in a prestressed elastic film on a Winkler foundation',
        generated_by='run_all.py',
        units_note='normalized (dimensionless) quantities; ell = sqrt(h Ebar/k), Ebar = E/(1-nu^2); '
                   'gamma = G Ebar/(sigma0^2 ell); kI = KI/(sigma0 sqrt(ell)); nu = 0.3 for kI values',
        definitions=dict(ell='sqrt(h Ebar / k)',
                         gamma_A='tanh(H/(2 ell))',
                         gamma_B='2 tanh(H/(2 ell)) - tanh(H/ell) = 2 tanh^3(H/2ell)/(1+tanh^2(H/2ell))',
                         kI_A='sqrt((1-nu^2) gamma_A)  (KII = 0 by symmetry)',
                         kI_B='sqrt((1-nu^2) gamma_B)  (KII = 0 by symmetry)'),
        sampling='ell/H in {0.5, 1, 1.5, 2}  =>  H/(2 ell) in {1, 0.5, 1/3, 0.25}',
        rows=rows,
        critical_stress=dict(
            formula='sigma0c = sqrt(Ebar Gc/(ell gamma))',
            small_H=dict(case_A='sigma0c ~ sqrt(2 Ebar Gc/H)', case_B='sigma0c ~ 2 ell sqrt(Ebar Gc/H^3)'),
            large_H=dict(both='sigma0c -> sqrt(Ebar Gc/ell) (isolated-crack threshold)')),
        checks=checks,
        missing_work=[])

    with open('results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(json.dumps({k: v for k, v in results.items() if k != 'rows'}, indent=2)[:800])
    print('rows:')
    for r in rows:
        print(r)
    allpass = all(c['passed'] for c in checks)
    print('ALL CHECKS PASS' if allpass else 'SOME CHECKS FAILED')
    print('wrote results.json')


if __name__ == '__main__':
    main()
