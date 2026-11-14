# The physics the agents missed (13 Nov 2026)

Teng, 13 Nov 2026, on the draft of knowledge levels: "no, these are not good. We cannot let agents search. what is the
physics the agents missed". The draft (`LEVELS_DRAFT.md`) named a paper or gave its formulas; it is withdrawn. This note
states the missed physics in mechanics terms. Organizer only: it contains the key's coefficients.

## In one sentence

Tension on a crack plane is a scalar problem with one compliance; shear is a vector problem with two, and a curved
front mixes them. The three runs that failed treated the two shear modes as copies of tension joined by a rotation.

## The two compliances of shear

Slip on a plane in an isotropic solid is resisted differently in two directions. Slip that varies along its own
direction is edge-like (plane strain) and stiffer by the factor 1/(1 − ν) than slip that varies across its own
direction, which is screw-like (anti-plane). In Fourier space on the crack plane, traction = −(μ/2)[|k| I +
(ν/(1 − ν)) k k/|k|] slip for the two in-plane components, and −(μ/(2(1 − ν)))|k| for opening. Opening has no such
direction: its factor is the same for every k.

The same two compliances appear at the front: for a given stress intensity the slip behind the front is
(1 − ν) K2 across the front and K3 along it (in units of 4/μ times the square root of r/2π).

Along a straight front the two never meet: mode 2 is purely edge-like, mode 3 purely screw-like. A curved front brings
them together in two ways.

## 1. Through the waviness: the two shear modes get different weights

A wavy front makes the fields vary along the front as well as across it, so mode 2 is no longer purely edge-like and
mode 3 no longer purely screw-like. The nonlocal effect of the front shape, which for tension is K1 = K1⁰(1 + L[a]),
therefore enters with a weight that is different for each shear mode and depends on Poisson's ratio:

    (2 − 3ν)/(2 − ν) = 1 − 2ν/(2 − ν) for mode 2,      (2 + ν)/(2 − ν) = 1 + 2ν/(2 − ν) for mode 3.

They are both 1 at ν = 0, where slip in any direction is resisted alike, and their mean is 1 for every ν. At ν = 0.3
they are 0.647 and 1.353.

- Missed by Kimi (weight 1 for all modes) and by DeepSeek (weight −1 for both shear modes, "shear modes unstable").
- Found by gpt-6.1-sol (recalled) and by GLM (computed from the crack-plane operator above).

## 2. Through the tilt: the coupling of the shear modes is not a rotation

Where the front is tilted by the slope a′, its normal and tangent turn, and part of each shear mode becomes the other.
If the stress intensity factors turned like the components of the traction ahead of the front, the couplings would be
(1, 1). If they turned like the components of the slip behind it, they would be (1/(1 − ν), 1 − ν). Neither is right,
because the traction and the slip have different compliances in the two directions. The published couplings,

    2/(2 − ν) into mode 2,      2(1 − ν)/(2 − ν) into mode 3,

are the harmonic means of those two pictures (an arithmetic fact of the formulas; reading it as a compromise between
turning the traction and turning the slip is our interpretation). At ν = 0.3: 1.176 and 0.824.

- Missed by DeepSeek (1 and 1: turning the traction) and half missed by Kimi (1 and 0.824).
- GLM had Poisson's ratio in the coupling with the right magnitude and the wrong way round (0.824 and 1.176).

## The principle that decides it, and that no failing run used

The energy release rate is the derivative of one elastic energy with respect to the position of the front. A small
tilt of the front cannot change it to first order, or equivalently its response to the front shape is a symmetric
operator. With G ∝ (1 − ν) K2² + K3² for the shear part, this requires

    b3 = (1 − ν) b2

for the two couplings. Turning the traction violates it by +ν, turning the slip by −ν; the published pair satisfies it.
The three failing runs violate it by +0.12 (Kimi), +0.30 (DeepSeek, exactly the turning of the traction) and +0.60
(GLM), which is why each of them reports a sine term in G for the cosine front of Part A. The task asks for that sine
coefficient, so the check was available to every run at no cost. Only gpt-6.1-sol made it.

## By run

| Run | Two weights (waviness) | Two couplings (tilt) | Energy unchanged by a tilt |
|---|---|---|---|
| gpt-6.1-sol | right | right | checked, holds |
| GLM-5.3 | right, computed | Poisson effect the wrong way round | violated, not checked |
| Kimi K3.1 Pro | shear treated as tension | one right, one a plain rotation | violated, seen and explained away |
| DeepSeek V4.1 Flash | shear treated as tension with the sign reversed | plain rotation | violated, not checked |

## What this says about knowledge to give

If knowledge is added to the task, it should be physics of this kind, stated without a source and without the
coefficients: that shear has two compliances which a curved front mixes, and that the energy release rate cannot
depend on the tilt of the front. Whether such statements are enough is not known; no run has been made with them, and
any text added to the task is for Teng to write or approve.

## What is established and what is not

- The weights and couplings are those of Gao and Rice (1986), Eq. 15, also in Lazarus (2011). The two weights were
  reproduced by GLM's own computation from the crack-plane operator; the two couplings were not reproduced by any
  computation independent of the paper.
- The relation b3 = (1 − ν) b2 follows from the energy argument alone.
- The statements about harmonic means and about the size of each run's violation are arithmetic, checked on
  13 Nov 2026. The physical reading of the harmonic mean is ours.
