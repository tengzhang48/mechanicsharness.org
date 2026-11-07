# REPRODUCE.md — independent checker

Work only in this folder. Do not edit `reference_submission/`.

## Commands

Independent deciding-quantity check (full-width variational strip, natural total-traction faces; does not call the author’s analytic functions):

```
python3 check_independent.py
```

Replay of the **unchanged** submitted code (copy first; the script writes `reference.json` next to itself):

```
mkdir -p author_replay
cp reference_submission/reproduce.py author_replay/reproduce.py
python3 author_replay/reproduce.py
```

Requirements: python3, numpy, scipy, sympy (author replay). No network. Runtime of each command is a few seconds.

## Artifacts

- `DERIVATION.md` — independent strip/energy derivation
- `check_independent.py` / `check_independent_output.txt` / `check_independent_results.json`
- `author_replay/reproduce_output.txt` / `author_replay/reference.json` (live replay)
- `ANSWER.md`, `results.json`, `PROVISIONAL_REPORT.md`, `STATUS.md`

## What the independent check does

P1 (and lumped-mass) FEM of I[v]=v(L)-v(0)+∫(v'^2+v^2)/2 on the full period. Stationarity implies v''=v and natural BCs v'(0)=v'(L)=-1. Then gamma_A=-I(H/ell), gamma_B=-2 I(H/ell)+I(2H/ell). Convergence, limits, and C1–C3 controls are reported in `check_independent_output.txt`.

Shared assumptions with the reference: 1D mode-I Winkler wake, PROBLEM.md energy, total traction, translating-tip bookkeeping. Not a 2D tip-field solve and not physical validation of a real film.
