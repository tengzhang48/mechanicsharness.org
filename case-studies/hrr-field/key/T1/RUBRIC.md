# Rubric — 100 points

Equivalent normalizations and numerical methods are welcome. Grade the requested quantities in the problem's convention; code reuse and named historical solutions earn no points by themselves.

| Item | Points | Full credit and explicit partial credit |
| --- | ---: | --- |
| Stress exponent | 15 | 10 for p within 10⁻⁵ of 0.133333333333; 5 instead if within 0.001. Another 5 for relating the nonzero contour flux to the radial energy exponent; 2 for a correct scaling argument that is incomplete. |
| Mechanics and boundary conditions | 25 | 10 for an equilibrated angular field and strain compatibility (5 if correctly set up but not solved); 10 for plane-strain closure and the correct effective-stress/strain factors (5 for correct closure with one factor error); 5 for symmetry, free faces, and an explicit amplitude normalization (1–4 for a correct subset). |
| Tensile coefficient H | 25 | 20 for H within 0.005 of 1.893569725213; 10 instead if within 0.05, or 5 for a sound numerical procedure with a larger calculation error. Another 5 for obtaining the amplitude from the stated J with consistent units and convention; 2 for recognizing this calibration but mishandling a constant. |
| Energy coefficient C | 20 | 15 for C within 0.002 of 0.457305783974; 8 instead if within 0.02, or 3 for a valid integral with a larger arithmetic error. Another 5 for volume/radial integration and the energy-per-thickness convention; 2 for correct angular energy but an incorrect radial or thickness factor. |
| Verification | 15 | 10 for a meaningful second method, radius-independent contour calculation, or analytic limiting check, with its numerical agreement reported; 5 for a relevant proposed check without results. Another 5 for demonstrated endpoint/equation residuals or normalization convergence; 2 for describing these without numerical evidence. |

Report beside the total whether all three numbers are at full credit (26 Oct: R31 scored 55 with H and C wrong; `key/check_P015.py` applies the bands). For each numerical row, choose one of its alternative score levels, then add any separately stated method points. Wrongly adding elasticity to the stipulated limiting law, using a plane-stress invariant, or using a table at a neighboring integer n loses the affected mechanics and amplitude points; it does not erase independently correct scaling or verification work. No full-credit finite-distance physical stress prediction is being requested.
