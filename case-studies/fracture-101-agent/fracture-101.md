# Fracture Mechanics 101 — Agent Edition

**Version:** 0.4.0 (draft of 16 November 2026, after two outside reviews; not yet used in a scored run)
**Audience:** AI agents that solve mechanics problems by analysis and computation
**Purpose:** physics that carries over from one crack problem to the next. It is not an answer key.

**How to read it.** Each entry has three parts: the **statement**, the **reason** for it, and its **use**: what it lets
you set up or test. Each entry holds under the conditions it names. The file gives standard relations and general
principles; it gives no solved task. If a task asks you to derive something that is stated here, the derivation is
still yours to give. Solve the problem that is stated to you, and use these entries to set it up and to test your
answer.

## A. Setting up a crack problem

### A1. The local frame of a crack front

- **Statement.** At each point of a smooth crack front there are three directions: the normal to the crack surface;
  the direction in that surface normal to the front, pointing into the uncracked material; and the tangent to the
  front. Take them as a right-handed set, with the normal to the crack surface pointing from the face called minus
  to the face called plus. The three modes are defined in this local frame. Mode I opens the faces. Mode II slides them across the
  front. Mode III slides them along the front.
- **Reason.** At distances from the front that are small compared with its radius of curvature and with every other
  length of the problem, the leading singular field at an interior point of the front has three independent parts:
  two of plane strain (modes I and II) and one of antiplane shear (mode III), in the plane normal to the front there.
- **Use.** When the front is not straight, the frame changes from point to point. A shear load whose direction is
  fixed in space is seen differently from points where the front has different directions; the local mode II and
  mode III factors are not obtained by resolving alone (entry D3). Decide first in which frame each quantity of the
  problem is defined. The direction normal to the front is the direction of advance
  only if the crack is held to its surface; otherwise the path is part of the problem (entry B2).

### A2. Jumps, areas, units and signs

- **Statement.** The opening or the slip of a crack is the jump of the displacement between its two faces, not the
  displacement of one face. Crack area is counted once, on one face. A stress intensity factor K has the unit
  stress × √length. An energy release rate G and a fracture energy Γ are energies per unit crack area.
- **Reason.** A rigid motion moves a face without opening the crack. The jump is twice the displacement of one face
  only when the problem is symmetric about the crack plane.
- **Use.** Label the two faces and write the jump as the value on the plus face minus the value on the minus face. Write down
  once the positive directions of the axes, of the tractions and of the jumps, and the sign convention of every
  transform you use, and keep them (entry B5).

### A3. Elastic constants near a front

- **Statement.** For an isotropic linear elastic solid with Young's modulus E and Poisson's ratio ν, the shear modulus
  is μ = E/[2(1 + ν)]. The modulus that enters the in-plane fracture relations is E′ = E/(1 − ν²) in plane strain and
  E′ = E in plane stress.
- **Reason.** At an interior point of a smooth front, the strain along the front is of lower order than the singular
  strains in the plane normal to the front. The leading in-plane field is therefore one of plane strain.
- **Use.** Use the plane strain value at an interior point of a three-dimensional front. Where the front meets a free
  surface or an interface, the near-tip field is that of the corner or of the interface; do not put plane stress in
  its place. Plane stress belongs to a thin sheet treated as two-dimensional, at distances from the tip that are
  large compared with its thickness.

## B. Crack-tip fields and energy

### B1. A stress intensity factor is one amplitude, and it can be read in two ways

- **Statement.** In the region where the singular field dominates (entry B6), in a homogeneous isotropic solid whose
  crack faces carry no singular load, the stress component of each mode on the crack plane ahead of the front is K/√(2πr) to leading order: the normal stress
  for mode I, the shear stress across the front for mode II, the shear stress along the front for mode III. Behind
  the front the displacement jump of each mode is, to leading order, (8K/E′)√(r/2π) for modes I and II and
  (4K/μ)√(r/2π) for mode III. In plane strain 8/E′ = 4(1 − ν)/μ.
- **Reason.** The stress ahead and the jump behind are two parts of one singular solution of the elasticity
  equations. The factor between them is the compliance of the mode. In plane strain, antiplane shear has 1/(1 − ν)
  times the compliance of opening and of in-plane shear.
- **Use.** From any solution, obtain K from the stress ahead and again from the jump behind, each by a limit toward
  the front. The two must agree (entry B5). In a formula of the type K = Yσ√(πa), the factor Y depends on the
  geometry, the boundary conditions and the distribution of the load.

### B2. The energy release rate

- **Statement.** G is the decrease of the potential energy per unit increase of crack area, at fixed loading. The
  potential energy is the elastic strain energy minus the work of the applied loads. Γ is the energy the material
  needs per unit crack area, with the same convention for the area. For a crack that grows slowly and cannot heal,
  with a resistance that does not depend on the rate: G ≤ Γ everywhere on the front, and the front advances only
  where G = Γ.
- **Reason.** Under fixed boundary displacements the loads do no work and the strain energy falls as the crack grows.
  Under fixed forces the strain energy rises and the potential energy falls by the same amount. The relation to the
  stress intensity factors follows from the work needed to close the crack again over a small advance: half the
  stress that acted ahead of the front, times the jump that appears behind it. With entry B1, for a homogeneous
  isotropic solid in plane strain and a crack that advances in its own plane,

      G = (K_I² + K_II²)/E′ + K_III²/(2μ).

  The three modes add because each stress does work only on its own jump.
- **Use.** Under mixed loading add the three terms; do not add the stress intensity factors. G does not say how the
  resistance Γ depends on the mix of modes, and it does not choose the path: say whether the crack is held to its
  plane and which Γ applies. A crack that is free to choose its path does not in general stay in its plane under
  shear. G = Γ gives an equilibrium state. Whether it is stable depends on how G − Γ changes with advance, and that
  depends on what is held fixed: the load, the displacement, or a loading device in between. When the body between the loading
  and the tip dissipates energy, keep apart the energy that flows into the tip and the energy supplied from far away.

### B3. One energy

- **Statement.** Take an elastic body under conservative loading, with the prescribed forces or displacements held
  fixed, a crack that advances within one fixed surface, and faces that do not rub. Its potential energy at equilibrium is one function of the shape of the crack. G at each point of the front
  is minus the derivative of that function with respect to the crack area added at that point.
- **Reason.** Nothing is dissipated in the body, so the energy depends on the present crack shape and loading, not on
  how the crack reached that shape. (If the crack leaves its surface, the energy depends on the whole surface it has
  made, not on the front alone.)
- **Use.** Second derivatives of one function do not depend on the order in which they are taken. Add a small crack
  area at a point A and record the change of G at a point B. Then add the same area at B and record the change of G
  at A. The two changes are equal. For a planar front written as x = a(z), the response of G at z to a change of a at
  z′ is the same when z and z′ are exchanged. If the front is described by a
  coordinate other than its arc length, equal areas mean a weight in this test. A computed response of G to the
  front shape that does not have this symmetry contains an error. Apply the test to one loading and one equilibrium state.

### B4. Superposition and weight functions

- **Statement.** At a fixed crack geometry, with linear material and linear boundary conditions, the fields and the
  stress intensity factors of two loadings add. G is quadratic in them and does not add. A weight function gives the stress intensity factor at a point of the
  front caused by a unit load on the crack faces.
- **Reason.** The elasticity equations and the conditions on a fixed crack are linear.
- **Use.** State whether the unit load of a weight function acts on one face or is a pair of opposite forces, and how
  it is normalized. A weight function belongs to its geometry and boundary conditions; using it elsewhere needs a
  derivation. Superposition is lost when the crack shapes differ, when the faces come into contact, or when a
  cohesive law is nonlinear. In three dimensions a shear load on the faces can produce both shear modes along the
  front (entry D2).

### B5. Two tests of a computed stress intensity factor, and signs

- **Two readings.** The reading from the stress ahead and the reading from the jump behind must agree (entry B1).
- **One energy.** The response of G to the front shape must have the symmetry of entry B3.
- **Overall balance.** The two tests above do not see an error in the overall amplitude. Check also that the stress
  carried by the uncracked material balances the applied load.
- **What a failure means.** A disagreement shows an error: a factor, a sign, the elastic constants, a frame, a reading
  taken outside the region where the singular field dominates, or a numerical error. Find it; do not explain it
  away.
- **Signs.** A term that changes sign when a coordinate direction is reversed has a sign that follows from the
  conventions of the calculation. The definitions of the axes, the faces and the positive loads fix the overall
  sense. The two tests fix the sign of one contribution relative to another. Passing them is necessary; it does not
  by itself prove that every sign is right.

### B6. Where the singular field holds, and the next term

- **Statement.** The K field describes the solution in a range of distances from the front: larger than the zone
  where the material fails or yields, smaller than the crack size, the specimen size, the radius of curvature of the
  front and the wavelength of any variation along it. The next term of the field is a stress that stays finite at
  the tip, such as a uniform stress parallel to the crack.
- **Reason.** The singular term is the leading term of an expansion in the distance from the front. How fast the
  next terms become small depends on the problem.
- **Use.** Before using K alone, check that such a range exists and that the finite stresses are small in it. Keep
  the next term when it acts on the zone at the tip, on the stability of the path, or on a second material nearby.

### B7. Open faces

- **Statement.** The relations above assume that the crack faces do not touch.
- **Use.** Check that the computed opening is nowhere negative. If the faces touch, the gap is zero there, the
  contact pressure is compressive, and friction needs its own law. A negative mode I factor is not an admissible
  open crack, and its square in entry B2 has no meaning. Sliding with friction dissipates energy; entries B3 and B4
  then do not apply.

## C. Beyond linear elastic fracture

### C1. Cohesive and bridging zones

- **Statement.** A traction–separation law gives the traction between the faces as a function of their jump. For a
  jump that grows steadily in one direction until the faces are fully separated, the area under the law is the work
  of separation per unit area. A strength σ₀ and a work of separation Γ define a length of the order E′Γ/σ₀², with
  the modulus of the mode in question.
- **Reason.** Material that is failing still carries load over a zone. The size of that zone is set by the balance
  between the energy the elastic body can deliver and the strength.
- **Use.** If the zone is short compared with the crack and the specimen, the body outside sees an effective
  resistance, which can still depend on the mix of modes, the rate and the history. While a bridged wake is still building up, the
  resistance rises as the crack grows; it becomes constant only once the zone is fully developed. If the zone is not
  short, solve
  the elastic body and the law together. The length E′Γ/σ₀² is an estimate; the zone length is that length times a
  number that depends on the law and the geometry. The fracture energy at the tip and the work done in a bridged wake
  are different quantities: say which one a result refers to. For a zone ahead of a traction-free crack, the stress
  is finite at the leading end of the zone, and the opening at its trailing end is not zero.

### C2. Plasticity, J and power-law fields

- **Statement.** If the plastic zone is small, there can be a range of distances outside it, and inside the other
  lengths of the problem, in which the field is the elastic K field. For a homogeneous elastic material with a
  stored-energy function, linear or not, uniform in the direction of advance, in a two-dimensional field with a
  straight traction-free crack and no body force, the integral J does not depend on the contour and equals the energy released per unit crack area.
- **Reason.** The contour independence of J rests on those hypotheses. A plastic material is not elastic; J can still
  be used as an approximation when the loading grows in proportion, without unloading. That is the case of a crack
  that does not grow; a growing crack unloads the material behind it.
- **Use.** For the usual power-law hardening model, in which the plastic strain grows as the stress to the power n,
  the stress near the tip of a crack that does not grow varies as r^(−1/(n+1)) and the strain as r^(−n/(n+1)), outside the zone of large strains
  at the blunted tip. The exponents follow from a scaling argument. The angular functions, the amplitude and any
  eigenvalue come only from the boundary-value problem. State the multiaxial form of the hardening law, what n
  means, and whether the field is plane strain or plane stress.

### C3. Cracks in a fluid-saturated porous solid

- **Statement.** In a linear poroelastic solid the unknowns are the total stress, the pore pressure and the fluid
  content, and fluid flows down the gradient of the pore pressure. Just after a sudden loading the fluid has not
  moved and the solid responds with its undrained constants. Pore pressure then diffuses over a distance of the order
  √(ct), with c the diffusivity. Consider a crack that does not move. At any fixed time t > 0, close enough to the
  tip, the pore pressure is bounded, and the relations between the singular stress, the displacement jump and the
  energy release rate are those of the drained solid, with its Poisson ratio and not the undrained one. This holds
  whether the crack faces are
  kept at a fixed pressure or are impermeable.
- **Reason.** The undrained pore pressure of an opening or in-plane shear loading is singular at the tip. It cannot
  persist: the flow it drives would dissipate an unbounded amount of energy. The time needed for pressure to even
  out over a distance r is of the order r²/c, which is shorter than
  any given t when r is small enough. The pore pressure at the tip is then bounded and the singular field carries
  the drained constants. With impermeable faces the fluid still moves through the solid around the tip, and the
  pressure at the tip need not be zero. In pure antiplane shear there is no volume change and no pore pressure to
  begin with.
- **Use.** At early times, while √(ct) is small compared with the crack and the specimen, the field farther than
  about √(ct) from the tip, and from any boundary where the undrained field does not meet the pressure or flux
  condition (an impermeable crack face can be one), is still close to the undrained one. A stress intensity factor
  defined by a limit toward the tip at fixed time belongs to the inner region: in general it differs from the value
  of the undrained body, and it can depend on the condition on the faces. Obtain it from the coupled problem or by joining the two
  regions; do not take it from a solution that is drained everywhere or undrained everywhere. The limit toward the
  tip and the limit toward zero time need not commute. Give the mechanical condition and the pressure or flux
  condition separately on every boundary, the crack faces included. A moving crack brings a second length, set by
  its speed.

### C4. Viscoelastic fracture

- **Statement.** A viscoelastic solid responds through its relaxation and creep functions, in general one pair for
  shear and one for volume change. If it has a finite short-time modulus and a nonzero long-time modulus, the first is
  the larger. Around a crack that grows steadily at speed V, material at distance r from the tip is loaded over a time
  of the order r/V.
- **Reason.** Where r/V is short compared with the relaxation times the material responds with its short-time
  modulus; where it is long, with its long-time modulus; in between it dissipates energy.
- **Use.** Compare r/V with the relaxation times over the range of distances that the problem really has, from the
  size of the failure zone to the size of the specimen. A short-time region and a long-time region exist only where
  that comparison allows them. The dependence on the speed comes from these comparisons: with no length at the tip
  and no length of the specimen, the ratio of the energy supplied from far away to the energy reaching the tip
  cannot depend on the speed. The energy supplied from far away pays for the fracture energy at the tip and for the
  dissipation in the body; keep the two apart. Creep compliance and relaxation modulus are different functions, and
  their characteristic times need not coincide. A finite specimen limits the region that can dissipate.

### C5. Films and diffusion along boundaries

- **Statement.** A stress that exists before the crack, such as the residual stress of a film, is part of the total
  traction on a crack face. On a free face the total traction is zero.
- **Reason.** The condition on a face is on the total traction, not on the change of traction.
- **Use.** If the problem is solved for the change of stress, the faces carry the prescribed traction minus the
  traction of the prior stress; on a free face this is the opposite of the prior traction. Reduce a film on a
  substrate to a membrane, a plate or a foundation model only when the stated problem gives that model or the
  conditions for it; a thick elastic substrate responds nonlocally, and each reduction changes the characteristic
  length and the energy release rate. When matter diffuses along a boundary and opens it, transport and elasticity are
  coupled; if the linearized problem has a slowest decaying mode, its decay time comes from the coupled problem with
  its boundary and conservation conditions. A time built from dimensional analysis gives the scale, not that number.

## D. Curved crack fronts in three dimensions

### D1. Moving one part of a front changes the loading of all of it

- **Statement.** Moving part of a planar front changes the boundary-value problem, and the stress intensity changes
  along the whole front. For a straight front under tension, perturbed in its own plane with a wavelength small
  compared with the other lengths, a part that is ahead of its neighbors is loaded less and a part that lags is
  loaded more.
- **Reason.** The crack faces are free of traction and the ligament is intact; changing where one ends and the other
  begins changes the elastic field everywhere, most strongly nearby.
- **Use.** The response of the stress intensity factors to the front shape is an integral over the front. A local
  rule, or a line tension with one coefficient, needs its own justification. Before concluding on stability, include
  how the loading of the straight front itself changes as it advances in the given specimen.

### D2. Shear on a crack plane is a vector problem with two compliances

- **Statement.** Consider slip on a plane in an infinite homogeneous isotropic solid, varying sinusoidally in one
  direction. Slip parallel to the direction in which it varies is resisted 1/(1 − ν) times as strongly as slip
  perpendicular to it (more strongly when ν > 0). If the geometry and the loading do not vary along a straight
  front, the first kind is mode II and the second is mode III, and they are uncoupled. When the fields also vary
  along the front, the two shear components can be coupled.
- **Reason.** Slip parallel to its direction of variation compresses and stretches the material on the two sides, as
  an edge dislocation does; slip perpendicular to it does not, as for a screw dislocation. Opening separates from
  shear when the problem is symmetric about the crack plane.
- **Use.** The two shear modes do not respond to the front shape as tension does, nor as each other. The response of
  tension cannot be copied into the shear modes. Solve the two shear components together.

### D3. A front that turns is not a rotation of the stress intensity factors

- **Statement.** Where the direction of the front changes, the mode II and mode III factors in the local frame change
  for two reasons: the local axes turn, and in general the elastic field changes.
- **Reason.** One could keep the shear stress amplitudes of the straight front fixed in space and resolve them on the
  new axes. One could equally keep the slip amplitudes fixed in space and resolve those. Because the two shear modes
  have different compliances (entry B1), the two procedures give different local factors. Neither follows from the
  local frame alone. (Turning the whole problem, the body, the crack and the loading together, is a different
  operation: it turns everything consistently and changes nothing.)
- **Use.** Do not obtain the local shear factors of a curved front by resolving the factors of the straight front on
  the new axes, and do not replace the curved front by a straight one turned through an angle. Take them from the
  solution for the actual front and loading, and apply both tests of entry B5.

### D4. A front held back by tough regions

- **Statement.** Where the fracture energy Γ varies along the front, an equilibrium front has G = Γ on the parts that
  advance and G ≤ Γ on the parts that are held (entry B2). This assumes that the front remains one curve and leaves no
  unbroken islands behind.
- **Reason.** G depends on the whole shape of the front (entry D1), so the toughness at a point does not fix the
  position of the front at that point.
- **Use.** The shape is an unknown of the problem: solve for it together with G. An equilibrium shape need not exist,
  need not be stable, and need not be the one reached under a given loading history.

## E. Limits of these entries

Each entry holds under the conditions it names. Formulas for small perturbations of a front need small perturbations;
the local frame and the local singular field do not. The isotropic formulas change at an interface between two
materials, in an anisotropic solid, at a corner of the front, and where the front meets a surface. In a solid whose
properties vary smoothly, the leading field keeps its form, with the properties at the tip. Large strains, plastic unloading, a failure zone that is not small, inertia and growth out of the
plane each need their own equations; the energy statements of entries B2 and B3 remain valid for elastic bodies. In
each case start from the governing equations and the boundary conditions of the problem at hand.

---

**Content boundary.** This file contains physics only: no task, no reference value, no grading key, no source. Any
coefficient, eigenvalue or profile that a task asks for has to come from your own solution of that task.
