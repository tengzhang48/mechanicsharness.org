# Knowledge entry: Hutchinson (1968), the power-law crack-tip field

Hutchinson J.W. (1968), Singular behaviour at the end of a tensile crack in a hardening
material. J. Mech. Phys. Solids 16, 13–31. Corpus id P015.

Topic: crack-tip asymptotics (`../../TOPIC.md`). Read the topic first: the governing equations,
how the exponent is fixed, J and its conventions, the scope table and the checks are there.
This entry holds what is specific to the paper. Prototype of the new format (log P240); its
content is reorganized from the verified and frozen v1 unit
(`corpus/knowledge_v1/P015_hutchinson_1968`, four correction rounds, Codex verification) and
adds no new scientific claim.

## 1. The paper's problem

A central tensile crack of half-length L in an infinite plate under remote tension σ∞ normal to
the crack; traction-free faces; mode I. Total deformation theory with Ramberg–Osgood hardening,
written non-dimensionally as ε = σ + α σⁿ (lengths by L, stresses by the yield stress σ_y,
strains by ε_y = σ_y/E). Plane stress and plane strain. The near-tip region is treated as a
semi-infinite crack.

The paper finds the dominant singular field as r → 0, where the plastic term dominates, and its
amplitude under small-scale yielding from Rice's J integral.

## 2. Results, with how each was established here

| Result | Status in this library |
|---|---|
| s = (2n+1)/(n+1) for the Airy function φ = K r^s φ̃(θ), i.e. stress exponent p = 2 − s = 1/(n+1), the same in plane stress and plane strain | derived (topic §3) and reproduced: an independent solver recovers all eight cases to 3.4×10⁻¹⁴ |
| Angular fields from the plane-stress and plane-strain eigenvalue ODEs (paper eqs. 13 and 31) with the face conditions | reproduced; the figures (3a/b, 5a/b) only qualitatively |
| Iₙ (eq. 24) and the energy coefficient Eₙ (eq. 34), Table 1 | reproduced: Table 1 within 0.27% (I) and 0.85% (Eₙ); see `reference_data.json` |
| Small-scale-yielding amplitude: K = (1/α)^{1/(n+1)} (π/I)^{1/(n+1)} (σ∞)^{2/(n+1)} (plane stress, eq. 27); the same with α replaced by α/(1 − ν²) (plane strain, eq. 32) | taken from J matching; consistent with the reproduced I |
| The tensile stress ahead of the tip is larger in plane strain than in plane stress, so plane-strain toughness is inferred to be **lower** (§7, eq. 50, Fig. 7) | taken; a field comparison, not itself a fracture criterion |
| Piecewise-linear material (eqs. 35–45): plane-stress amplitude factor √(E_tan/E) (eq. 42) | taken, not computed; the author calls the single-term version "suspicious" |

## 3. Conventions (state them whenever you use a number from here)

- Normalization: **max σ̃ₑ = 1** over θ. Iₙ and Eₙ in Table 1 depend on it.
- Exponent notation: the paper uses s for the Airy function, φ ∝ r^s; stresses go as r^(s−2).
- K here is the HRR amplitude (load intensity); I is a dimensionless angular normalization
  integral. They are different objects. The plane-strain factor (1 − ν²) enters only through
  the remote elastic J = π(1 − ν²)σ∞² (non-dimensional), not the relation between J and K.
- Small-scale yielding fixes only the amplitude: J from an outer elastic contour is πσ∞²
  (plane stress) or π(1 − ν²)σ∞² (plane strain). The local exponents do not need it.

## 4. Reference data

`reference_data.json`: Table 1 as printed (I, (π/I)^{1/(n+1)}, Eₙ for n = 3, 5, 9, 13, both
constraints), with transcription notes, and independently computed values to 12 digits (Codex,
verification round 3). Use them to check your own solver in the paper's case; they are not
answers to other problems.

## 5. The paper as printed: misprints and transcription notes

1. Eq. (22), traction-work line: printed exponent r^((n−1)(s−2)) and no α; it must be
   α K^(n+1) r^((n+1)(s−2)) (dimensions, and consistency with eqs. 23–24). Eqs. (23)–(27) use
   the correct form.
2. Eq. (37): the hydrostatic term is printed ((1 − 2ν)/6) σ_KK; the complementary energy needs
   σ_KK².
3. Eq. (44): the angular differential is printed dσ instead of dθ.
4. Table 1 in the text layer (OCR) is unreliable: plane-stress I(13) is 2.87 on the page (OCR:
   2.37); plane-strain Eₙ(9) is 1.48 on the page (OCR: 1.43). Read the page image.
5. The printed eq. (24) bracket, including the −sin θ(…) term, is correct. A "first-principles"
   J giving I ≈ 5.65 (plane stress, n = 3) came from dropping polar-basis derivatives
   (topic §4), not from the paper.

## 6. Leaving the scope (specific to this paper)

The topic's table applies. In addition:

- Deformation theory is justified only for monotonic, nearly proportional loading. The paper
  examines this only afterwards (§8: stress proportions differ considerably between the
  elastic and plastic zones on the same ray) and calls for a check; it is not established.
- The dominant term describes only a small region inside the plastic zone. Small-scale
  yielding (Rice 1966 suggested applied stress below about 1/3, at most 1/2, of yield) is needed
  for the amplitude formula, not for the local exponents.
- n = 1 is linear elasticity, where the plastic/elastic split of the paper degenerates;
  E_tan = 0 (perfect plasticity) is a singular limit not described by this solution.
- Smooth-yielding materials (e.g. aluminium); not the Lüders-band behaviour of mild steel (§7).

## 7. Reference code (optional)

`reference_code/hrr.py` is the v1 solver, unchanged (SHA-256 in `reference_code/SHA256SUMS`,
identical to the frozen v1 copy, whose 70 tests pass there). Its output passes every topic check
(`../../checks/test_checks.py`). It is **a crack-only reference**: the wedge angle is
fixed at π, mode I, pure power law. Functions starting with `_` are internal. For any other
geometry, loading or law, write your own solver and test it with the topic's checks; use this
code, if at all, to compare with in the paper's own case.
