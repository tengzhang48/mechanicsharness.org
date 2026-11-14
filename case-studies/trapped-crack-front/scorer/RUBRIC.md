# Physics rubric of the trapped-front benchmark

Approved by Teng on 13 Nov 2026, after the first four runs: "answer close is one thing. missing a physics is another.
If one agent misses a physics, they need to lose points"; "We just need to update the rubric"; and, on the rubric
below, "please document. I think it is ok."

## The rule

`score.py` decides whether a number is right, by the rule stated in the task: within 1% of the reference, or within
0.005 if that is larger. That rule and the 76 checks are unchanged. `rubric.py` then gives points by the physics the
numbers test, 100 in all.

**The response, 70 points.** Seven pieces of physics, 10 points each. Each is tested by one or two numbers that isolate
it, and is earned only if all of them pass.

| Piece of physics | Tested by | Points |
|---|---|---|
| Tension: a part of the front that lags is loaded more | C1, S1 | 10 |
| Energy release rate in mixed mode and the load level of equilibrium | K10_scaled | 10 |
| Shear across the front: its own weight, which depends on Poisson's ratio | C2 | 10 |
| Shear along the front: its own weight, which depends on Poisson's ratio | C3 | 10 |
| The slope of the front turns shear along it into shear across it | S2 | 10 |
| The slope of the front turns shear across it into shear along it | S3 | 10 |
| The energy release rate of a wavy front, with no term in the slope | CG, SG | 10 |

**The outcome, 30 points.** The equilibrium front in the tough bands (front_peak_to_peak and the sixteen values of f):
10 points. The stress intensity factors along it (the three sixteen-point tables, kappa3_max, s_kappa3_max): 20 points.
The outcome needs every piece of the response, so it is counted only in a run whose response is complete; then each
part earns its points in proportion to its numbers that pass.

## What follows from it

- A run that misses a piece of physics loses that piece's points and the 30 points of the outcome that rests on it.
- A number that is close but outside the tolerance earns nothing. Closeness can be reported as a description of an
  answer; it is not credit.
- With the whole response right, a slip in one table value costs a fraction of a point (0.4 for one of the fifty
  numbers of the stress intensity tables), not the outcome.
- The count of right numbers (of 76) is still reported beside the points. The two can order runs differently: a run
  can have more right numbers and less of the physics.

## How it was tested

`check_controls.py`, part 5: the reference earns 100 points, also when rounded to three decimals; each of seven wrong
theories, solved consistently, loses exactly the pieces it misses and the whole outcome (10 to 50 points); one wrong
table value in an otherwise right answer gives 99.6.

## The runs of 13 Nov 2026 under it

| Run | Points | Numbers right | Pieces of the response missing |
|---|---|---|---|
| codex_1 (gpt-6.1-sol) | 100 | 76 of 76 | none |
| glm_1 (GLM-5.3) | 40 | 8 of 76 | both slope couplings; the energy release rate of the wavy front |
| kimi_1 (Kimi K3.1 Pro) | 30 | 12 of 76 | both shear weights; the coupling into shear across the front; the energy release rate of the wavy front |
| deepseek_1 (DeepSeek V4.1 Flash) | 20 | 3 of 76 | everything but tension and the load level |

Per run: `../runs/<run>/score/rubric.md`.

## What it does not change

The rubric was written after these runs: the agents were told only the pass rule. The task text is unchanged. Stating
the rubric in the task would be a new version of the task text (not done).
