"""Assemble an agent's workspace for this benchmark: the task text only, then check that nothing else is in it.
Organizer only.

usage: python3 make_workspace.py WORKSPACE_DIR             the task text as it is (the runs of 13 Nov 2026)
       python3 make_workspace.py WORKSPACE_DIR reminder    the same text with the reminder of task/REMINDER.md
(the folder must not exist)

The reminder (Teng, 13 Nov 2026, after the first runs: "I think we need to just remind the AI agents to carefully check
the shear modes on curved crack front") is one paragraph, placed after the opening paragraph of the task and before
section 1. Nothing else differs: the script checks that removing that paragraph gives back the task text byte for byte.
In both cases the workspace holds one file, PUBLIC_TASK.md.
"""
import os, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__))
GIVEN = ['task/PUBLIC_TASK.md']
MARK = '## 1. The body and the crack\n'          # the reminder goes in front of this heading
ws = sys.argv[1]
variant = sys.argv[2] if len(sys.argv) > 2 else 'none'
if variant not in ('none', 'reminder'):
    sys.exit('unknown variant %r (none or reminder)' % variant)
if os.path.exists(ws):
    sys.exit('%s exists; refusing to overwrite' % ws)
os.makedirs(ws)
for rel in GIVEN:
    shutil.copy2(os.path.join(HERE, rel), os.path.join(ws, os.path.basename(rel)))
if variant == 'reminder':
    base = open(os.path.join(HERE, 'task/PUBLIC_TASK.md'), encoding='utf-8').read()
    reminder = open(os.path.join(HERE, 'task/REMINDER.md'), encoding='utf-8').read().strip()
    assert base.count(MARK) == 1 and reminder and '\n' not in reminder
    block = reminder + '\n\n'
    text = base.replace(MARK, block + MARK)
    assert text.replace(block, '', 1) == base and text.count(reminder) == 1
    with open(os.path.join(ws, 'PUBLIC_TASK.md'), 'w', encoding='utf-8') as fh:
        fh.write(text)
found = sorted(os.path.relpath(os.path.join(d, f), ws) for d, _, fs in os.walk(ws) for f in fs)
expected = sorted(os.path.basename(r) for r in GIVEN)
if found != expected:
    sys.exit('workspace content %s differs from %s' % (found, expected))
print('workspace', ws, 'with', ', '.join(found), '(variant: %s)' % variant)
