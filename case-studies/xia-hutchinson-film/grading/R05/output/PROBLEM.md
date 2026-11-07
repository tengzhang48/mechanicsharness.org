# Energy release for periodic cracks in an elastic film

Draft for Teng's review. This is the proposed source-faithful reproduction
task; source solutions and the proposed geometry extension are excluded.

## Mathematical model

An infinite isotropic elastic film of uniform thickness h is attached to a
local elastic foundation. Its Young's modulus is E, Poisson ratio is nu,
and its reference state has uniform equibiaxial tensile stress sigma0 I.
Assume h, E, k and sigma0 are positive and -1 < nu < 1/2. Use nu = 0.3
for the numerical values requested below. The response is quasistatic and
uses the plane-stress constitutive law. The substrate remains elastic,
bonding is maintained, and all cracks go through the film thickness.

Let u(x1, x2) be the in-plane displacement relative to this reference state,
and let Greek indices range over 1 and 2. Define

\[
\epsilon_{\alpha\beta}
=\tfrac12(u_{\alpha,\beta}+u_{\beta,\alpha}),\qquad
\Delta\sigma_{\alpha\beta}
=\frac{E}{1-\nu^2}
\left[(1-\nu)\epsilon_{\alpha\beta}
+\nu\epsilon_{\gamma\gamma}\delta_{\alpha\beta}\right],
\qquad
\sigma_{\alpha\beta}=\sigma_0\delta_{\alpha\beta}
+\Delta\sigma_{\alpha\beta}.
\]

The foundation exerts force per unit film area f_alpha = -k u_alpha,
where k has units Pa/m. Equilibrium and the crack-face condition are

\[
h\sigma_{\alpha\beta,\beta}-ku_\alpha=0,
\qquad \sigma_{\alpha\beta}n_\beta=0.
\]

Displacement may jump across each crack. Define the change of total
potential energy in a film region Omega, relative to its prestressed state,
by

\[
\Delta\Pi
=\int_\Omega
\left[
h\sigma_0\epsilon_{\alpha\alpha}
+\frac h2\Delta\sigma_{\alpha\beta}\epsilon_{\alpha\beta}
+\frac k2 u_\alpha u_\alpha
\right]dA.
\]

The reference prestress and model parameters are fixed during an increment
of crack advance. For bookkeeping, Ebar denotes E/(1 - nu^2). Energy
release G is energy per advancing crack area: divide release per unit
in-plane advance per tip by h. It has units J/m^2. Use conventional local
plane-stress stress intensity factors; in this model
G = (KI^2 + KII^2)/E.

## Two crack geometries

Case A consists of semi-infinite cracks on x2 = mH, x1 < 0, for every
integer m, with H > 0. Each tip is at x1 = 0. The displacement increment
tends to zero far ahead, and the relaxed cracked wake is reached far
behind. The field is periodic in x2 with period H. All tips advance an
equal distance in the positive x1 direction.

Case B already contains infinite cracks on x2 = 2mH for all x1. Additional
semi-infinite cracks lie on x2 = (2m + 1)H, x1 < 0, with tips at x1 = 0.
Far ahead, the film is in the equilibrium state of the existing infinite
cracks. Far behind, the additional cracks split that wake. The field is
periodic in x2 with period 2H. Only the additional tips advance, equally
in the positive x1 direction; the old cracks remain fixed.

Every face of every crack is traction free. For each case the far-ahead
and far-behind states are the bounded equilibria of the stated model.

## Requested solution

1. Determine the decay length ell of the crack-induced wake field and
   derive the relevant displacement and stress fields far ahead and far
   behind the tips. State the conditions satisfied at every strip edge.
2. Determine GA(H) and GB(H), per advancing tip. Justify the energy
   accounting, including the period width, moving-tip count and thickness
   normalization. Determine which fracture modes the geometries permit and
   report KI/(sigma0 sqrt(ell)).
3. For a positive film toughness Gc, give the critical reference stress for
   advance in each case. Derive the H/ell -> 0 and H/ell -> infinity limits
   at fixed h, E, nu and k, and explain their mechanical meaning.
4. Report GA Ebar/(sigma0^2 ell), GB Ebar/(sigma0^2 ell), and both
   KI/(sigma0 sqrt(ell)) at ell/H = 0.5, 1, 1.5 and 2. These are normalized
   sampling locations, not prescribed dimensional material data.
5. Give reproducible checks of the deciding quantities and explain what
   errors each could detect. State which physical conclusions are supported
   by this reduced model and what requires additional information or
   validation.

## Deliverables

Provide ANSWER.md with the derivation, assumptions, results and unresolved
questions; results.json with clearly named normalized quantities and their
sampling locations; and any code with a single reproduction command.
Report checks with actual outputs and declared tolerances. Identify shared
assumptions between verification routes.

If a necessary specification is unclear, record the question explicitly;
do not silently change the problem. A provisional report should state
completed quantities, current checks and unresolved issues, with partial
work retained.
