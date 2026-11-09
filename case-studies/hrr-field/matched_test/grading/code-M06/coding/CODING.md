- **B1 existence: present_right; before.** The finite-J argument selects the admissible exponent s = 2/9, consistent with the reference; the later face residual corroborates the assumed angular field. Evidence: `submission/PROCESS.md:5`, "J finiteness forces s=1/(n+1)=2/9 (stress exponent p=1/(n+1))"; `submission/PROCESS.md:7`, "got f'(pi)=1e-16 (face condition automatic)".

- **B2 limiting cases: present_right; after.** The same angular solver and flux calculation reproduce the known n = 1 field and normalization, matching the reference. Evidence: `submission/PROCESS.md:8`, "linear-limit check n=1 reproduced f=2 sin(theta/2), f(pi)=2, I_J=pi analytically and numerically".

- **B3 applicability: present_right; after.** The final answer identifies the small-strain intermediate regime and its inner limit, and excludes outer geometry and a fracture criterion; these are physically appropriate limits for the stationary leading field. Evidence: `submission/ANSWER.md:81`, "The field is a small-strain intermediate asymptotic; very near the tip it becomes self-inconsistent (finite strains / process zone)"; `submission/ANSWER.md:81`, "the outer geometry and the actual fracture criterion are outside this calculation."

- **B4 degenerate cases: absent; unclear.** Neither file examines coinciding exponents, a threshold where the singularity appears or disappears, or proximity to another change of structure. There is no evidence quote or assignable timing.

Timing uses the first angular ODE integration (`submission/PROCESS.md:7`) as the start of the main numerical work; the n=1 check follows it but precedes the amplitude and requested outputs. The applicability discussion is documented only in the final answer, and the n=1 comparison does not examine a structural degeneracy.
