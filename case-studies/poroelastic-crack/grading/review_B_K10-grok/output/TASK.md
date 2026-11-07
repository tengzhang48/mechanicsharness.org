# Task: settle a contested result

Problem `PROBLEM.md` (a suddenly pressurized crack in a Biot poroelastic solid) was solved by agents.
The key (`KEY.md`) follows the source paper, Detournay & Cheng 1991 (`paper/`; text layer and PDF):
for loading 1, K(0+) = [(1 - nu_u)/(1 - nu)] sqrt(pi L), i.e. 0.75 sqrt(pi L) for nu = 0.2, nu_u = 0.4
(their Eq. 29, which they call a "heuristic argument", p. 1650; compare their numerical Fig. 9).

Three independent runs of another model (`gpt_answers/MP-37, MP-38, MP-39`: ANSWER.md and code)
instead obtain K(0+) ≈ 0.807-0.808 sqrt(pi L), defined as lim(t->0+) lim(r->0+) sqrt(2 pi r) sigma_yy,
arguing that a drained tip layer sits inside the undrained outer field and that a modulus-ratio
estimate is not enough (a Wiener-Hopf factorization; an energy identity with a pressure-gradient term).
Two of them also dispute the paper's small-time volume coefficient (2.44 vs 1.354); the third agrees
with the paper (1.354).

Your job: decide, with evidence someone else can check, which is right for the stated problem.
1. Derive K(0+) independently (any valid route: matched asymptotics, Wiener-Hopf, a careful numerical
   solution resolving the tip layer, or a known limit such as the high-speed limit of a moving crack in
   a poroelastic solid if you can justify that it applies). State the order of limits you use.
2. Check the GPT derivation and code step by step; run their code on a copy in ./scratch/ if useful.
3. Assess whether the paper's numerical method (six quadratic elements over 2L, Stehfest inversion,
   p. 1655) can resolve the tip layer at the smallest times of Fig. 9.
4. Do the same for the small-time volume coefficient (1.354 vs 2.44).

Write `VERDICT.md`: for each of the two quantities, your value with uncertainty, which side it supports
(paper / GPT / neither), the decisive argument, and a check that could have failed. Be explicit about
anything you could not settle. Keep `PROCESS.md` (times pasted from `date -u +%H:%M:%S`). Work only in
this folder; no web search or fetch. Time limit: 90 minutes. Never stop or kill a process you did not
start; kill only by recorded PID; never use pkill, killall or pgrep-based loops.
