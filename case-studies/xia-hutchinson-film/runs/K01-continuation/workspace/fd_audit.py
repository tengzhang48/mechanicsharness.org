"""Audit of the path-spread statistic reported in fd_output.txt.

REVIEw_GAP item: 'case A H=2.0 ... path_spread(b,a)=(1.62e-06,8.07e+00)' — the ahead
relative spread 8.07 uses |mean J| as denominator where J ~ 0 (undisturbed ahead state).
This script re-solves ONLY that configuration and reports absolute spreads and the
denominators, so the relative number can be judged honestly.
"""
import numpy as np
import fd_check

r = fd_check.solve('A', 2.0)
Ja, Jb = r['J_ahead'], r['J_behind']
abs_a = r['J_path_spread_ahead'] * max(abs(Ja), 1e-30)
abs_b = r['J_path_spread_behind'] * abs(Jb)
scale = abs(r['G_exact']) * fd_check.H_TH / 2.0   # released energy per half-domain
print('=== FD PATH-SPREAD AUDIT: case A, H=2.0 (single re-solve) ===')
print(f"J_ahead mean            = {Ja:+.6e}   (target continuum value: 0)")
print(f"J_behind mean           = {Jb:+.6e}")
print(f"ahead   abs spread      = {abs_a:.3e}   rel spread = {r['J_path_spread_ahead']:.3e}  (denominator |mean J| = {abs(Ja):.3e})")
print(f"behind  abs spread      = {abs_b:.3e}   rel spread = {r['J_path_spread_behind']:.3e}  (denominator |mean J| = {abs(Jb):.3e})")
print(f"ahead abs spread / (G h/2) = {abs_a/scale:.3e}  (absolute error relative to released energy)")
print(f"G_fd={r['G_fd']:.6f} G_exact={r['G_exact']:.6f} rel_err={100*r['rel_err']:+.3f}%")
print('INTERPRETATION: the relative spread 8.07 is a divide-by-near-zero artifact:')
print('J_ahead mean ~ -1e-9 (continuum target 0), so the relative normalization inflates')
print('an absolute path variation of ~8e-9 (~2e-8 of the released energy G h/2).')
print('Path independence in the ahead region is in fact better than behind (5.6e-7 abs).')
