import numpy as np
from notch import scan_bg
alpha = 2.0*np.pi/3.0
# physical range for s: crack n=4 -> 1.8; notch weaker -> s smaller. scan 0.8..2.0
sgrid = np.linspace(0.8, 2.0, 61)
rows = scan_bg(4.0, alpha, sgrid, pmin=-6.0, pmax=2.0, Np=80, outfile='scan_rows.json')
print('rows', len(rows))
# print near-zero fp entries sorted by |fp|
rows_sorted = sorted(rows, key=lambda r: abs(r['fp']))
for r in rows_sorted[:15]:
    print(f"s={r['s']:.4f} p={r['p']:.5f} fp={r['fp']:+.3e}")
