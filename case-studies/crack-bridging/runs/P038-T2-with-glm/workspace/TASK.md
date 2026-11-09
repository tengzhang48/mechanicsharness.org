# Task

Solve the problem below. Write your complete answer to `ANSWER.md` in this folder: the requested
numbers or decisions, your method, and your checks. You may write and run code in this folder
(`python3` with numpy and scipy). Budget: 45 minutes; write `ANSWER.md` early, update it, and make
sure it is complete by 40 minutes. Keep `PROCESS.md`: one line per step, each starting with a time
pasted from `date -u +%H:%M:%S` (never type a time).

Rules: use only this folder (including `./knowledge/`) and your own knowledge. Do not open files
outside this folder, and do not use web search or web fetch. The session log will be checked.

Provided material: `./knowledge/` contains a knowledge unit about a mechanics paper (an entry with a
"Using this unit" section, an extraction, a report, Python code and its tests). Use it as you see fit.

---

# Separation between unequal strips

Two long parallel strips, each $b=25.0\ \mathrm{mm}$ wide, are joined along a
flat weak plane. Their longitudinal elastic moduli and thicknesses are:

| Strip | Modulus | Thickness |
|---|---:|---:|
| Upper | $110\ \mathrm{GPa}$ | $0.60\ \mathrm{mm}$ |
| Lower | $190\ \mathrm{GPa}$ | $0.90\ \mathrm{mm}$ |

A $20.0\ \mathrm{mm}$ starter cut extends from the left end along the joining
plane. The bonded part far ahead is long and remains straight. Equal and opposite
opening couples of magnitude $M$ are applied at the two free ends, with **no end
forces**. $M$ denotes the total moment on **one strip**, not a moment per unit
width or the sum of the two moments.

The strips remain linear elastic with small rotations. Axial stretching and
transverse shear deformation are negligible; lateral contraction is free, so use
the quoted longitudinal moduli for bending. As separation grows beyond the starter
cut, the new separated faces retain connectors. They exert only closing normal
traction, with monotonic opening data

\[
t(\delta)=\sigma_b\left(1-\frac{\delta}{\delta_c}\right),
\quad 0<\delta<\delta_c,\qquad
\sigma_b=0.20\ \mathrm{MPa},\quad \delta_c=0.20\ \mathrm{mm}.
\]

A connector becomes permanently inactive upon reaching $\delta_c$. The gap
$\delta$ is the full relative opening. Advancing the separation front consumes
$\Gamma_0=1.40\ \mathrm{J/m^2}$ in addition to connector work, per area of the
joining plane, counting the pair of new faces once. Treat the intrinsic front
process as spatially negligible. There is no other damage or tangential face force.

Increase the couples quasistatically from initiation. Let $L$ be the distance
from the advancing front to the original starter-cut root, and let $x$ run from
that front ($x=0$) toward the root ($x=L$). Determine, when a **newly retained
connector first reaches $\delta_c$ and becomes inactive**:

1. $M$ in N·m on each strip;
2. $L$ in mm;
3. $\delta(L/2)$ in µm.

Briefly justify your mechanical description and boundary conditions, and give one
check. Each number must be within **±0.5%** of the reference. Budget: 45 minutes;
Python with NumPy and SciPy is available.
