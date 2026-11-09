# Builder notes

Source: G. Bao and Z. Suo (1992), *Remarks on Crack-Bridging Concepts*,
Applied Mechanics Reviews 45, 355–366, corpus P038. Only the supplied
[text layer](paper.txt) and local [page images](pages/) were used; no web search
or fetch was performed. Formula signs and factors were checked against the
scans on printed pages 359, 362 and 364.

Each probe stands alone. Supply only its TASK_PROBLEM.md to a solver;
REFERENCE.md, RUBRIC.md, reference.py and reference_results.json are grader
materials. No statement identifies the source, cites an equation number, tells
the solver which numerical method to choose, or relies on another probe.
Every statement supplies physical data and the opening/area conventions.

## What each probe tests

| Probe | Transfer being tested | Relation to the source | Expected effort |
|---|---|---|---|
| T1, parameter transfer | A finite bridged flaw at new dimensional inputs, retaining intrinsic front energy; simultaneous tip-growth and connector-failure conditions | Section 3.2, finite-crack relations on printed p. 362 (equations 45–46), with the nonzero front-work criterion from the general framework | Two scalar conditions reducible to one positive root, or an equivalent elasticity computation; within 45 minutes |
| T2, variant | Adapt splitting mechanics to dissimilar arm stiffnesses, a softening law and nonzero front work; compute geometry-dependent zone/profile | Section 3.4, splitting-beam energy balance and geometric zone-length effects on p. 364, combined with the softening law of Section 2.3 | An unknown-length fourth-order boundary problem, or an exact solution; within 45 minutes using NumPy/SciPy |
| T3, applicability | Distinguish localized elastic bridging, failure of bulk elasticity, and opening-constrained partial use of a valid local law | The small/large-scale distinction of Sections 1–3, the local-law/component coupling principle, and the splitting-versus-tunneling discussion of Section 3.4 | Three short physical judgments; no full tunnel or plasticity simulation |

T1's numerical inputs and result are not tabulated or stated in the source.
Its intrinsic front work is positive, so using the zero-front-energy strength
curve alone is insufficient. The retained-zone length is 1.70 times the original
half-flaw: replacing it with a material toughness and the original traction-free
crack overpredicts the requested stress by 23.67%.

T2 changes both material arrangement and local law from the equal-arm,
constant-traction splitting illustration. For equal-width arms, the opening
compliance is the sum of the two separate bending compliances, not the
compliance of an averaged strip. The triangular curve area fixes the limiting
couple through energy balance but cannot fix the zone length or midpoint
opening. Intrinsic front work imposes nonzero tip curvature. A pure end couple
is transmitted unchanged through the starter cut, so its 20 mm length does
not affect these three requested zone quantities.

T2 explicitly prescribes elastic strips with negligible transverse shear and
axial stretching, free lateral contraction and small rotations. Its reference
is the corresponding strip-bending idealization, not a finite-thickness
two-dimensional crack-tip correction. The computed retained length is 12 times
the thicker arm and the rear relative rotation is 0.0434 rad. These assumptions
make the problem definite and manageable; no governing equation, effective
stiffness, tip condition or characteristic root is supplied to the solver.

T3 does not announce that there is one case of each label. A permits the local
model and a small-zone reduction for developed resistance. B deliberately
violates the elastic-bulk premise by observed ligament-wide yielding; rejecting
the supplied complete description is appropriate, but a local relation can
still be retained inside an elastic-plastic extension. C is borderline in the
useful sense that a constrained elastic/local-law calculation remains possible,
but the full opening-test energy is unavailable. Its cap permits at most
0.726% of the full curve work, while that partial-work scale is still 25.8%
of the intrinsic front work. Neither a large full-energy toughening claim nor
an automatic zero-fiber-effect conclusion is supported.

## Reference checks and reproducibility

The programs [T1/reference.py](T1/reference.py),
[T2/reference.py](T2/reference.py) and [T3/reference.py](T3/reference.py)
are self-contained. They require only Python, NumPy and SciPy, with no
source-paper access.
From this folder:

    python T1/reference.py
    python T2/reference.py
    python T3/reference.py

Captured outputs are stored beside each program. Assertions implement the
reported numerical checks:

- T1: analytic finite-crack reduction versus independently integrated
  stress-intensity/opening Green functions, checked again by the uniform-pressure
  elliptical opening and zero-front-work/large-crack limits.
- T2: unknown-length boundary-value solver versus an exact
  hyperbolic/trigonometric solution and energy balance, plus constant-traction
  and zero-front-work special cases.
- T3: adaptive curve quadratures versus analytic antiderivatives, with energy
  length, yielding and capped-opening limits checked separately. The qualitative
  judgments are supported by physical premise/limit checks, not presented as
  experimentally validated predictions.

Observed relative agreement is below $10^{-13}$ for T1/T2 main outputs and
below $10^{-15}$ for T3 work quadratures. Requested tolerances (0.3% and 0.5%)
are deliberately much wider than numerical uncertainty. T3 has no mandatory
numerical output. No reference quantity is left unchecked.

## Grading guidance

Each rubric totals 100 points and states partial-credit rules. Score the
mechanical reasoning and outputs; source recall and code reuse are not required.
Accept equivalent finite-element, integral or exact approaches when they respect
the stated assumptions and conventions.

Full-gap versus one-face displacement, half-flaw versus total slit length,
front work versus connector work, and total moment versus moment per width are
explicit in the statements. Apply numerical penalties to propagated numerical
errors while preserving sound model/adaptation credit.

For T3, “no to full-energy transfer, yes to a constrained local-law calculation”
is the intended conditional answer for C, regardless of its headline label.
Qualitative reasoning can earn full credit. The partial area under C's curve
must not be called an automatically dissipated toughness or an actual tunnel
energy release rate. Lack of connector failure does not imply lack of connector
force. Blanket refusal and unjustified unconditional application have equal
case-score caps; supported adaptations and requests for truly missing inputs
are rewarded.

The chronological [PROCESS.md](PROCESS.md) uses times obtained directly from
the required UTC date command. Temporary inspection/validation work stays in
[scratch/](scratch/).
