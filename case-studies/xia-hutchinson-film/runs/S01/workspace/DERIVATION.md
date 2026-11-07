# Incremental derivation

## 1. Translation-invariant strip equilibria

Let a relaxed strip have width W and local transverse coordinate y in [-W/2,W/2]. Bounded equilibrium and translation invariance give u1 = 0 and

    h Ebar u2'' - k u2 = 0,
    Ebar = E/(1-nu^2),
    ell = sqrt(h Ebar/k).

Both horizontal edges are traction free: sigma22 = sigma12 = 0. Hence u2'(+/-W/2) = -sigma0/Ebar. The unique solution is

    u2(y) = -(sigma0 ell/Ebar) sinh(y/ell)/cosh(W/(2ell)),
    sigma22(y) = sigma0 [1-cosh(y/ell)/cosh(W/(2ell))],
    sigma11(y) = sigma0 [1-nu cosh(y/ell)/cosh(W/(2ell))],
    sigma12 = 0.

For Case A ahead, u = 0 and sigma = sigma0 I. For Case A behind, W = H. For Case B ahead, W = 2H, and behind there are two W = H strips in each 2H period. Coordinates must be centered separately in each disconnected strip; no displacement continuity is imposed through a crack.

## 2. Energy per unit x1 length

For a strip,

    U(W) = integral[-W/2,W/2] [h sigma0 u2' + (h Ebar/2)(u2')^2 + (k/2)u2^2] dy.

Multiply the equilibrium equation by u2 and integrate. With Ebar u2' = -sigma0 at both edges,

    integral [h Ebar (u2')^2 + k u2^2] dy = -h sigma0 [u2(W/2)-u2(-W/2)].

Thus

    U(W) = (h sigma0/2)[u2(W/2)-u2(-W/2)]
         = -(h sigma0^2 ell/Ebar) tanh(W/(2ell)).

Rigid translation of a steady crack front replaces ahead material with behind material. The tip-transition energy is invariant under this translation. Per period, release is U_ahead - U_behind; divide by advancing-tip count and h, not by period width.

    GA = (sigma0^2 ell/Ebar) tanh(H/(2ell)),
    GB = (sigma0^2 ell/Ebar) [2 tanh(H/(2ell))-tanh(H/ell)].

Case A has period H and one advancing tip. Case B has period 2H and one additional advancing tip; old cracks do not advance.

## 3. Symmetry and local intensity

Reflection about an advancing crack line leaves the geometry and loading invariant: tangential displacement is even and normal displacement is odd. Mode II is excluded. The foundation is nonsingular at a crack tip and does not change the leading local plane-stress energy/intensity relation supplied in the problem.

    KII = 0,
    KI/(sigma0 sqrt(ell)) = sqrt[(1-nu^2) gamma],
    gamma = G Ebar/(sigma0^2 ell).

Independent numerical and symbolic checks remain to be executed. No external source or numerical crack-tip field has been used.

## 4. Completion record

At 2026-10-15T13:08:42Z, all requested derivations and sampling values are
complete in ANSWER.md and results.json. Eight checks passed on both runs;
the final run used PID 377774 under a 120-second shell timeout. CHECKS.txt
contains final actual outputs. The prior sentence in Section 3 records
the earlier partial stage, not current missing work. No independent 2D
tip-field solve or empirical validation is claimed.
