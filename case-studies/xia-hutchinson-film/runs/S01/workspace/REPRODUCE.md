# Reproduction

Run from this folder:

```bash
timeout 120s python check.py
```

The command requires local NumPy, SciPy and SymPy, reads no scientific input files or external data, and writes CHECKS.txt and results.json incrementally. All dimensionless sampling values use nu = 0.3 and the specified ell/H values. The script prints its PID and UTC time. The checks solve only the far-field strip model, not a separate two-dimensional crack-tip problem.
