"""Assemble an agent's workspace for this benchmark: the task text only, then check that nothing else is in it.
Organizer only. usage: python3 make_workspace.py WORKSPACE_DIR [VARIANT]   (the folder must not exist)
VARIANT: none (default) or no_definitions: task/VARIANT_NO_DEFINITIONS.md, given as PUBLIC_TASK.md."""
import os, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__))
VARIANTS = {'none': 'task/PUBLIC_TASK.md', 'no_definitions': 'task/VARIANT_NO_DEFINITIONS.md'}
variant = sys.argv[2] if len(sys.argv) > 2 else 'none'
if variant not in VARIANTS:
    sys.exit('unknown variant %r; known: %s' % (variant, ', '.join(VARIANTS)))
GIVEN = [VARIANTS[variant]]
ws = sys.argv[1]
if os.path.exists(ws):
    sys.exit('%s exists; refusing to overwrite' % ws)
os.makedirs(ws)
for rel in GIVEN:
    text = open(os.path.join(HERE, rel), encoding='utf-8').read()
    assert '<!--' not in text, 'organizer notes in %s' % rel
    shutil.copy2(os.path.join(HERE, rel), os.path.join(ws, 'PUBLIC_TASK.md'))
found = sorted(os.path.relpath(os.path.join(d, f), ws) for d, _, fs in os.walk(ws) for f in fs)
expected = ['PUBLIC_TASK.md']
if found != expected:
    sys.exit('workspace content %s differs from %s' % (found, expected))
print('workspace', ws, 'with', ', '.join(found))
