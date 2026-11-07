"""Build the HTML summary of the multiphysics round (log P390) from summary_data.json (runs, automatic
checks, blind grades) and the Phi(eta) values of B/DIMENSIONAL_ANALYSIS.md. All numbers on the page
come from those files; text statements cite the documents they summarize.
Usage: python3 grade_v2.py profile; python3 summary_data.py; python3 make_summary_html.py
Writes SUMMARY.html here and a copy in <documents>."""
import json, os, html, shutil
from collections import Counter, defaultdict

H = os.path.dirname(os.path.abspath(__file__))
D = [r for r in json.load(open(f"{H}/summary_data.json")) if r["problem"] in ("A", "B", "D", "E")]   # this round only; gel and battery runs are reported elsewhere
AG = [("gptaz", "GPT-6.1-sol", "opencode, Azure"), ("k3az", "Kimi K3", "Cline, Azure"),
      ("glm", "GLM-5.3", "opencode"), ("dsazT", "DeepSeek-V4.1-Flash", "Cline, Azure, reasoning on (rerun)"),
      ("dsaz", "DeepSeek-V4.1-Flash", "reasoning OFF (configuration error; superseded)")]
PR = [("A", "Tip of a fluid-driven crack", "Detournay 2016; Garagash & Detournay", 50),
      ("B", "Suddenly pressurized crack, poroelastic", "Detournay & Cheng 1991", 50),
      ("D", "Cracks in piezoelectric ceramics", "Suo, Kuo, Barnett & Willis 1992", 50),
      ("E", "Grain-boundary diffusion wedge in a thin film", "Gao, Zhang, Nix, Thompson & Arzt 1999", 170)]
esc = html.escape
by = defaultdict(list)
for r in D:
    by[(r["agent"], r["problem"])].append(r)

def cell(agent, prob):
    runs = sorted(by.get((agent, prob), []), key=lambda r: r["repeat"] or 0)
    if not runs:
        return '<td class="notrun">not run</td>'
    rows = []
    for r in runs:
        if r["parts"]:
            sq = "".join(f'<span class="sq {v}" title="({k}) {v}"></span>' for k, v in sorted(r["parts"].items()))
            kc = f'{r["key_checks"][0]}/{r["key_checks"][1]}' if r["key_checks"] else "–"
            meta = f'<span class="meta">key {kc} · {r["minutes"]} min · {r["errors"]["decisive"]} decisive</span>'
        else:
            sq = '<span class="sq none"></span>' * 3
            meta = '<span class="meta">no answer</span>'
        rows.append(f'<div class="run"><span class="rl">{esc(r["label"])}</span><span class="sqs">{sq}</span>{meta}</div>')
    return '<td>' + "".join(rows) + "</td>"

totals = {}
for a, *_ in AG:
    c = Counter()
    for r in D:
        if r["agent"] == a:
            c.update(r["parts"].values() if r["parts"] else ["none"] * 3)
    totals[a] = c

# Phi(eta) chart: exact values (two independent factorizations; B/DIMENSIONAL_ANALYSIS.md) against candidates
PHI = [(0.02, 0.985267), (0.1, 0.931222), (0.2, 0.872917), (1 / 3, 0.807593), (0.4, 0.779170), (0.52, 0.733709),
       (0.6, 0.706831), (0.8, 0.649112), (1.0, 0.601850)]
W_, H_, ml, mr, mt, mb = 640, 360, 58, 150, 18, 46
X = lambda e: ml + e * (W_ - ml - mr)
Y = lambda v: mt + (1.0 - v) / 0.55 * (H_ - mt - mb)
def curve(f):
    return " ".join(f"{X(i / 100):.1f},{Y(f(i / 100)):.1f}" for i in range(101))
svg = [f'<svg viewBox="0 0 {W_} {H_}" role="img" aria-label="Phi of eta">']
for v in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
    svg.append(f'<line class="grid" x1="{ml}" x2="{W_ - mr}" y1="{Y(v):.1f}" y2="{Y(v):.1f}"/>'
               f'<text class="tick" x="{ml - 8}" y="{Y(v) + 4:.1f}" text-anchor="end">{v:.1f}</text>')
for e in (0, 0.2, 0.4, 0.6, 0.8, 1.0):
    svg.append(f'<text class="tick" x="{X(e):.1f}" y="{H_ - mb + 18}" text-anchor="middle">{e:g}</text>')
svg.append(f'<text class="axl" x="{(ml + W_ - mr) / 2:.0f}" y="{H_ - 6}" text-anchor="middle">η = (ν_u − ν)/(1 − ν_u)</text>')
svg.append(f'<text class="axl" x="14" y="{(mt + H_ - mb) / 2:.0f}" text-anchor="middle" transform="rotate(-90 14 {(mt + H_ - mb) / 2:.0f})">K_tip / K_outer at t = 0⁺</text>')
lines = [("c1", lambda e: 1 / (1 + e), "1/(1+η): D&amp;C 1991, Hui 2013"),
         ("c2", lambda e: (1 + e) ** -0.5, "(1+η)^−½: energy matching"),
         ("c3", lambda e: (1 + e) ** -0.735, "(1+η)^−0.735: Yu et al.")]
for cls, f, lab in lines:
    svg.append(f'<polyline class="{cls}" fill="none" points="{curve(f)}"/>'
               f'<text class="lab {cls}t" x="{X(1.0) + 6:.1f}" y="{Y(f(1.0)) + 4:.1f}">{lab.split(":")[0]}</text>')
for e, v in PHI:
    svg.append(f'<circle class="pt" cx="{X(e):.1f}" cy="{Y(v):.1f}" r="4"/>')
svg.append(f'<text class="lab ptt" x="{X(1.0) + 6:.1f}" y="{Y(0.60185) + 18:.1f}">exact (ours)</text>')
svg.append(f'<line class="mark" x1="{X(1/3):.1f}" x2="{X(1/3):.1f}" y1="{mt}" y2="{H_ - mb}"/>'
           f'<text class="tick" x="{X(1/3) + 4:.1f}" y="{mt + 12}">problem B</text>')
svg.append("</svg>")

challenges = [
    ("B (a), (c)", "K(0⁺) depends on the order of limits; the face-only volume growth misses a tip term", "GPT ×2", "ok", "Verified (Codex, Grok, our FE); key corrected"),
    ("B (b)", "No decaying 2D steady field: take the long-time limit pointwise", "GPT, GLM, K3", "ok", "Right; it is how the paper's limit is meant"),
    ("B (b)", "…so a uniform p means free swelling and K₂(∞) = V₂(∞) = 0", "K3 ×3, GLM", "bad", "Fails (ours): the pressurized region stays confined; interior stress −η·p at any radius"),
    ("A (b)", "ℓ is defined only up to a numerical factor", "all three agents", "amb", "Statement ambiguity; define ℓ next time"),
    ("A (c)", "The exact opening at 10 m differs from the viscosity asymptote by 5–20%", "K3 ×2, DeepSeek ×2", "bad", "Not supported (ours): the next term of the expansion changes it by 9 × 10⁻⁵"),
    ("D (b)", "Charge control reverses the electric part of 𝒢", "GLM, K3, DeepSeek", "open", "Likely fails (ours): 𝒢 is the same under charge or potential control; not computed"),
    ("E (c)", "The relaxation criterion is unspecified", "GLM", "amb", "Statement ambiguity"),
    ("E", "Initial traction vs free surface at t = 0; kernel sign; π normalization", "K3, DeepSeek", "open", "Open"),
]
cls_lab = {"ok": "verified / right", "bad": "fails", "amb": "ambiguity", "open": "open"}

def tr(cells):
    return "<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"

ans = sum(r["answered"] for r in D)
page = f"""<title>Multiphysics Crack Round</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Source+Serif+4:opsz,wght@8..60,600;8..60,700&display=swap">
<style>
:root {{ --bg:#f5f4ef; --surface:#fffefb; --ink:#1c2226; --muted:#5b646c; --rule:#dcd8cd; --accent:#1e5d78;
  --right:#2f7d50; --partly:#c28a1d; --wrong:#b23b3b; --contested:#6c51a3; --none:#a3a8ad; --soft:#ebe8df; color-scheme:light }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#13171a; --surface:#1a1f23; --ink:#e7e4dc; --muted:#9ba4ac;
  --rule:#2d3439; --accent:#73b6d2; --right:#5fc289; --partly:#e2b552; --wrong:#e47070; --contested:#a990e2; --none:#626a71; --soft:#222a2f; color-scheme:dark }} }}
:root[data-theme="dark"] {{ --bg:#13171a; --surface:#1a1f23; --ink:#e7e4dc; --muted:#9ba4ac; --rule:#2d3439; --accent:#73b6d2;
  --right:#5fc289; --partly:#e2b552; --wrong:#e47070; --contested:#a990e2; --none:#626a71; --soft:#222a2f; color-scheme:dark }}
body {{ background:var(--bg); color:var(--ink); font:15px/1.6 "IBM Plex Sans", system-ui, sans-serif; }}
.wrap {{ max-width:1080px; margin:0 auto; padding-inline:clamp(16px,4vw,40px); padding-block:40px 64px; display:grid; gap:44px }}
h1,h2,h3 {{ font-family:"Source Serif 4", Georgia, serif; text-wrap:balance; margin:0 }}
h1 {{ font-size:clamp(28px,4vw,40px); line-height:1.15 }}
h2 {{ font-size:23px; margin-bottom:10px }} h3 {{ font-size:17px }}
p {{ margin:0; max-width:68ch }} .lede {{ color:var(--muted); max-width:72ch }}
.kicker {{ font:500 12px/1 "IBM Plex Mono", monospace; letter-spacing:.08em; text-transform:uppercase; color:var(--accent) }}
section {{ display:grid; gap:14px; min-width:0 }}
.facts {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(230px,1fr)); gap:14px }}
.fact {{ background:var(--surface); border:1px solid var(--rule); border-radius:8px; padding:14px 16px; display:grid; gap:6px }}
.fact b {{ font-family:"Source Serif 4", Georgia, serif; font-size:17px }}
.fact span {{ color:var(--muted); font-size:14px }}
.scroll {{ overflow-x:auto; border:1px solid var(--rule); border-radius:8px; background:var(--surface) }}
table {{ border-collapse:collapse; width:100%; font-size:14px }}
th,td {{ text-align:left; vertical-align:top; padding:10px 12px; border-bottom:1px solid var(--rule) }}
th {{ font-weight:600; background:var(--soft) }} tr:last-child td {{ border-bottom:0 }}
.map td {{ min-width:190px }} .map th:first-child {{ min-width:150px }}
.agent b {{ display:block }} .agent span {{ color:var(--muted); font-size:12.5px }}
.run {{ display:grid; grid-template-columns:44px auto; column-gap:8px; align-items:center; margin:3px 0 }}
.rl {{ font:12px "IBM Plex Mono", monospace; color:var(--muted) }}
.sqs {{ display:inline-flex; gap:3px }}
.sq {{ width:16px; height:16px; border-radius:3px; display:inline-block }}
.sq.right {{ background:var(--right) }} .sq.partly {{ background:var(--partly) }} .sq.wrong {{ background:var(--wrong) }}
.sq.contested {{ background:var(--contested) }} .sq.none {{ background:transparent; border:1.5px dashed var(--none) }}
.meta {{ grid-column:2; font:11.5px "IBM Plex Mono", monospace; color:var(--muted); font-variant-numeric:tabular-nums }}
.notrun {{ color:var(--muted); font-style:italic }}
.legend {{ display:flex; flex-wrap:wrap; gap:14px; font-size:13px; color:var(--muted) }}
.legend span {{ display:inline-flex; align-items:center; gap:6px }}
.num {{ font-family:"IBM Plex Mono", monospace; font-variant-numeric:tabular-nums }}
.two {{ display:grid; grid-template-columns:minmax(0,1.1fr) minmax(0,1fr); gap:24px; align-items:start }}
@media (max-width:820px) {{ .two {{ grid-template-columns:1fr }} }}
figure {{ margin:0; background:var(--surface); border:1px solid var(--rule); border-radius:8px; padding:12px; min-width:0 }}
figure svg {{ width:100%; height:auto; display:block }}
figcaption {{ font-size:13px; color:var(--muted); margin-top:8px }}
svg .grid {{ stroke:var(--rule); stroke-width:1 }} svg .tick, svg .axl {{ fill:var(--muted); font:12px "IBM Plex Sans", sans-serif }}
svg .c1 {{ stroke:var(--wrong); stroke-width:2; stroke-dasharray:6 4 }} svg .c2 {{ stroke:var(--partly); stroke-width:2; stroke-dasharray:2 3 }}
svg .c3 {{ stroke:var(--accent); stroke-width:1.6 }} svg .pt {{ fill:var(--ink) }} svg .mark {{ stroke:var(--muted); stroke-dasharray:2 4 }}
svg .lab {{ font:12px "IBM Plex Sans", sans-serif }} svg .c1t {{ fill:var(--wrong) }} svg .c2t {{ fill:var(--partly) }} svg .c3t {{ fill:var(--accent) }} svg .ptt {{ fill:var(--ink) }}
.tag {{ font:500 11.5px "IBM Plex Mono", monospace; padding:2px 7px; border-radius:99px; white-space:nowrap; border:1px solid currentColor }}
.tag.ok {{ color:var(--right) }} .tag.bad {{ color:var(--wrong) }} .tag.amb {{ color:var(--partly) }} .tag.open {{ color:var(--muted) }}
ol, ul {{ margin:0; padding-left:20px; display:grid; gap:6px; max-width:72ch }}
code {{ font:13px "IBM Plex Mono", monospace; background:var(--soft); padding:1px 5px; border-radius:4px; overflow-wrap:anywhere }}
.foot {{ color:var(--muted); font-size:13px; display:grid; gap:6px }}
</style>
<div class="wrap">
<header style="display:grid;gap:12px">
  <div class="kicker">crack-front pilot · multiphysics round · 2026-10-08</div>
  <h1>Multiphysics Crack Round</h1>
  <p class="lede">Four coupled-field crack problems (fluid-driven, poroelastic, piezoelectric, diffusive) were given closed-book
  to four AI agents, three runs each (GPT on B only; DeepSeek twice, first with reasoning off by mistake): {len(D)} runs, {ans} with an answer. Each answer was rerun, checked against a
  derived key, and graded blind against a profile, not a single score. All numbers below come from
  <code>summary_data.json</code> (the run logs, automatic checks and blind grades).</p>
</header>

<section>
  <div class="kicker">At a glance</div>
  <div class="facts">
    <div class="fact"><b>Standard methods: mostly right</b><span>On A, D and E, K3 and GLM usually reached the key (K3: all
    D parts right in 3 of 3 runs and all A key checks passed in 3 of 3; GLM all A key checks in 2 of 3). DeepSeek, rerun with reasoning on, joined them on A and D.</span></div>
    <div class="fact"><b>Knowing the condition: only GPT</b><span>On B, GPT found the drained crack tip in 3 of 3 runs.
    None of the 9 K3, GLM and DeepSeek runs did: 8 took the elastic value K(0⁺) = 1 (one in other units), one gave 1/(1 − ν).</span></div>
    <div class="fact"><b>A challenge that held</b><span>GPT's challenge to the source paper was confirmed by Codex, Grok and our FE.
    The key was corrected (0.75 → 0.8076). The literature itself is split on this value.</span></div>
    <div class="fact"><b>Harness: three corrections</b><span>{ans} of {len(D)} runs answered. DeepSeek ran with reasoning off
    (Cline without <code>--thinking</code>), so its results understate it. GLM's three missing or late runs hit the provider's 5-hour
    usage limit. The other two failures were Cline stops (compaction; a text-only turn).</span></div>
  </div>
</section>

<section>
  <h2>The map: every run, every part</h2>
  <p class="lede">One square per question part, (a) (b) (c), as graded blind. Below each run: automatic key checks passed, minutes used
  (budget 50 min; E 170 min) and the grader's count of decisive errors.</p>
  <div class="legend"><span><i class="sq right"></i>right</span><span><i class="sq partly"></i>partly</span><span><i class="sq wrong"></i>wrong</span>
  <span><i class="sq contested"></i>contested (a checkable challenge, graded both ways)</span><span><i class="sq none"></i>no answer</span></div>
  <div class="scroll"><table class="map">
    <tr><th>Agent</th>{"".join(f'<th>{p} · {esc(t)}<br><span style="font-weight:400;color:var(--muted);font-size:12.5px">{esc(s)}</span></th>' for p, t, s, _ in PR)}</tr>
    {"".join(f'<tr><td class="agent"><b>{esc(n)}</b><span>{esc(h)}</span><br><span class="num">right {totals[a]["right"]} · partly {totals[a]["partly"]} · wrong {totals[a]["wrong"]} · contested {totals[a]["contested"]}</span></td>' + "".join(cell(a, p) for p, *_ in PR) + "</tr>" for a, n, h in AG)}
  </table></div>
</section>

<section>
  <h2>Problem B: the drained crack tip</h2>
  <div class="two">
    <div style="display:grid;gap:12px">
      <p>A crack in a fluid-saturated solid is loaded suddenly. In the bulk the response is undrained. At the tip the undrained
      pore pressure is singular, so the local drainage time r²/c goes to zero: at every t &gt; 0 a zone of size √(ct) has drained.
      K is defined at r → 0, so K(0⁺) = Φ·K_outer, not K_outer.</p>
      <p><b>Why the others missed it</b> (from the run files): a formulation with no place for the tip (undrained field plus a face layer);
      a standard theorem that made K = 1 look exact; checks that verified the elastic solution but could not fail on the assumption;
      early stops (K3 11–22 min, DeepSeek 6–16 min of 50). GLM spent 90k–150k reasoning tokens on the same chain.</p>
      <p><b>What GPT did:</b> wrote the exact transform formulation, in which the drained and undrained limits are the two ends of one
      function; started from the definition of K; sent adversarial checks to subagents; refined before concluding.</p>
      <p><b>The dimensional analysis behind it:</b> the global group ct/L² is not enough; near the tip the field needs r/√(ct), and the two
      orders of limits are its two ends. The inner problem is self-similar, so K_tip/K_outer = Φ(η) does not depend on time.</p>
    </div>
    <figure>{"".join(svg)}
      <figcaption>The tip factor Φ for drained faces against the coupling η. Dots: exact values from two independent factorizations
      (GPT's and Codex's codes). Lines: the published approximations. The paper's 0.75 for B is the dashed line at the marked η.</figcaption>
    </figure>
  </div>
</section>

<section>
  <h2>Challenges the agents raised</h2>
  <p class="lede">The grader records every checkable challenge to the problem, the key or a published result. "Ours" means Opus's assessment;
  only B's was verified by independent computation.</p>
  <div class="scroll"><table>
    <tr><th>Part</th><th>Challenge</th><th>Raised by</th><th>Status</th></tr>
    {"".join(tr([esc(p), esc(c), esc(b), f'<span class="tag {k}">{cls_lab[k]}</span> {esc(s)}']) for p, c, b, k, s in challenges)}
  </table></div>
</section>

<section>
  <h2>What this suggests, and what comes next</h2>
  <ol>
    <li>On coupled-field problems with a standard method, agents do well, as in milestone 1 for elastic problems.</li>
    <li>The hard part is the condition: where a simplification stops holding. On B the non-GPT runs made the same error twice,
    a limit taken over the wrong length: the undrained tip in (a), and free swelling in (b) (K₂(∞) = 0 in 6 of 9 runs).</li>
    <li>A challenge is the first step to an insight, but it has to be checked: B's challenge held; A's and the free-swelling one
    did not (ours).</li>
    <li><b>The entry test (done, automatic checks; blind grades deferred):</b> on the gel problems G1 and G2, DeepSeek found the drained
    tip in 3 of 3 runs with a short dimensional-analysis entry and in 0 of 3 without, on both problems. Two easier problems of the same
    family (M, C), which point at the thin region directly, were solved without help.</li>
  </ol>
</section>

<footer class="foot">
  <p><b>Caveats.</b> DeepSeek ran without reasoning in this round (<code>HARNESS_NOTES.md</code>); a rerun with reasoning on is
  needed before comparing it. Three runs per agent and problem; GPT ran only B. The blind grader is Codex with gpt-6.1-sol (the same model family
  as the GPT solver); the planned second-family grading of a sample is not done yet. The run-time model labels are kept in a
  coordinator-only file.</p>
  <p><b>Sources.</b> <code>multiphysics_round/ROUND_RECORD.md</code>, <code>B/WHY_DRAINED_TIP_MISSED.md</code>,
  <code>B/DIMENSIONAL_ANALYSIS.md</code>, <code>B/CASE_STUDY_GPT_CHALLENGE.md</code>, <code>B/FE_CHECK.md</code>,
  <code>gel/LITERATURE_CHECK.md</code>, <code>../MILESTONE_2/MILESTONE_2_PLAN.md</code>; built by <code>make_summary_html.py</code>.</p>
</footer>
</div>
"""
open(f"{H}/SUMMARY.html", "w").write(page)
shutil.copy(f"{H}/SUMMARY.html", os.path.expanduser("<documents>/Multiphysics_Round_Summary_2026-10-08.html"))
print("written", len(page), "bytes")
