# Fracture Mechanics 101, agent edition (version 0.4.0)

`fracture-101.md` is a general text of fracture physics written for AI agents. It is the companion of
[Fracture Mechanics 101](https://mechanicsharness.org/case-studies/fracture-101.html), which is for human readers. Each entry gives a statement, the reason
for it and its use: what it lets an agent set up or test. It contains no solved task, no answer key and no names of
source papers.

**How it was written.** ChatGPT wrote a first draft on 15 November 2026. Claude rewrote it on 16 November 2026 and
revised it after reviews by two other AI models (GPT-6.1 Sol and Gemini). Both reviews found no wrong formula and
asked for corrections, which were made. Teng Zhang directed the work. No expert in fracture mechanics has checked
the text line by line.

**What has been tested: one run.** On the [trapped crack front](https://mechanicsharness.org/case-studies/trapped-crack-front.html) case, Kimi K3.1 Pro,
working with a DeepSeek reviewer inside the run, was given this file and got all 76 requested numbers right. The same
setup without the file got 9 of 76.

**What that does not show.**

- The run was not repeated. On the same day the same agent scored between 20 and 40 points of 100 in seven other
  runs of this task, each with a different addition to the task text.
- The text was written after agents had failed on this task, and two entries (B1 and D3) state the physics behind
  errors made there. They are general statements. Still, the run shows a repair on this task; it does not show that
  the text helps on a problem it was not written from.
- A test on a second problem, a kinked crack, was begun. The agent solved it without the file, so that problem
  cannot show what the file adds.

**If you use the file in a test,** give it unchanged (SHA-256
`a64b83d4d795c71fc2eeaa5e9af7124855cc8e87c3d87d80638d7de01a7e5778`), record everything else the agent received, and
compare with the same setup without the file.

The header of the file says "not yet used in a scored run". That line was written before the run and is left as it
is, so that the file is byte for byte what the agent received.
