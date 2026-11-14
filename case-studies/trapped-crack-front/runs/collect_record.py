"""Copy the record of one finished run into runs/NAME/ of this folder (the part of a run that is kept in git).

Organizer only. Copies from RUNS_ROOT/NAME: the run's facts (benchmark commit, model, prompt, versions, workspace
hashes, task variant if any), AUDIT.md, score/ (score.py and rubric.py output), the grade (RUNS_ROOT/grading/PACKET.grade,
if a packet name is given; as grade/ if an attempt was accepted, as grade_not_accepted/ if none was), the submission (every file of the workspace except the task text, tmp/
and caches), and messages.md: the agent's visible messages and tool steps, cut from the event stream (opencode.jsonl,
which stays outside git). The stream holds no private reasoning.

A Codex run has codex.log instead of opencode.jsonl: its messages are cut by audit_codex_log.py, and the final message
that Codex wrote into the workspace (codex_final.md) is kept as final_message.md, not as part of the submission.
Another benchmark keeps its records in its own folder: the environment variable RECORD_DIR names it (default: this folder).

usage: python3 -I collect_record.py RUNS_ROOT NAME [PACKET]
"""
import json, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEST = os.path.abspath(os.environ.get('RECORD_DIR') or HERE)
FACTS = ['benchmark_commit.txt', 'model.txt', 'prompt.txt', 'opencode_version.txt', 'codex_version.txt',
         'workspace_hashes.txt', 'task_variant.txt', 'AUDIT.md']
SKIP = {'PUBLIC_TASK.md', 'tmp', '__pycache__', 'codex_final.md'}
MAX_FILE = 2_000_000     # a larger file of a submission is listed, not copied


def short(x, n=160):
    x = ' '.join(str(x).split())
    return x if len(x) <= n else x[:n] + ' ...'


def messages(stream, name):
    out = ['# Visible messages and tool steps of %s, in order' % name, '',
           'Cut from the event stream of the run (opencode.jsonl, kept outside git). The stream holds the agent\'s messages and its',
           'tool calls; it holds no private reasoning. Long tool inputs are shortened here.', '']
    n_msg = 0
    for line in open(stream, encoding='utf-8', errors='replace'):
        try:
            j = json.loads(line)
        except Exception:
            continue
        part = j.get('part') or {}
        if j.get('type') == 'text' and part.get('text', '').strip():
            n_msg += 1
            out += ['', '## Message %d' % n_msg, '', part['text'].strip(), '']
        elif j.get('type') == 'tool_use':
            st = part.get('state') or {}
            inp = st.get('input') or {}
            first = next(iter(inp.values()), '') if isinstance(inp, dict) else inp
            if isinstance(inp, dict) and inp.get('prompt'):        # a call to a sub-agent: what it was asked matters
                first = short(inp['prompt'], 700)
            out.append('- `%s` (%s): `%s`' % (part.get('tool'), st.get('status'), short(first, 700).replace('`', "'")))
            if part.get('tool') == 'subagent' and st.get('output'):   # and the start of what the sub-agent answered
                out.append('  - answer (start): `%s`' % short(st['output'], 900).replace('`', "'"))
    return '\n'.join(out) + '\n', n_msg


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    root, name = sys.argv[1], sys.argv[2]
    packet = sys.argv[3] if len(sys.argv) == 4 else None
    src, dst = os.path.join(root, name), os.path.join(DEST, name)
    if os.path.exists(dst):
        sys.exit('%s exists; refusing to overwrite' % dst)
    if not os.path.isfile(os.path.join(src, 'score', 'score.json')):
        sys.exit('no score in %s: score the run first' % src)
    os.makedirs(os.path.join(dst, 'submission'))
    for f in FACTS:
        if os.path.isfile(os.path.join(src, f)):
            shutil.copy2(os.path.join(src, f), os.path.join(dst, f))
    shutil.copytree(os.path.join(src, 'score'), os.path.join(dst, 'score'))
    notes = []
    ws = os.path.join(src, 'workspace')
    for d, dirs, files in os.walk(ws):
        dirs[:] = sorted(x for x in dirs if x not in SKIP)
        for f in sorted(files):
            rel = os.path.relpath(os.path.join(d, f), ws)
            if rel in SKIP:
                continue
            size = os.path.getsize(os.path.join(d, f))
            if size > MAX_FILE:
                notes.append('%s (%d bytes): not copied, larger than %d bytes' % (rel, size, MAX_FILE))
                continue
            os.makedirs(os.path.dirname(os.path.join(dst, 'submission', rel)), exist_ok=True)
            shutil.copy2(os.path.join(d, f), os.path.join(dst, 'submission', rel))
    if notes:
        open(os.path.join(dst, 'submission_not_copied.txt'), 'w').write('\n'.join(notes) + '\n')
    stream = os.path.join(src, 'opencode.jsonl')
    n_msg = None
    if os.path.isfile(stream):
        text, n_msg = messages(stream, name)
        open(os.path.join(dst, 'messages.md'), 'w', encoding='utf-8').write(text)
    elif os.path.isfile(os.path.join(src, 'codex.log')):
        text = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'audit_codex_log.py'), src, 'messages'],
                              capture_output=True, text=True, check=True).stdout
        n_msg = text.count('\n## Message ')
        open(os.path.join(dst, 'messages.md'), 'w', encoding='utf-8').write(text)
    if os.path.isfile(os.path.join(ws, 'codex_final.md')):
        shutil.copy2(os.path.join(ws, 'codex_final.md'), os.path.join(dst, 'final_message.md'))
    graded = False
    if packet:
        g = os.path.join(root, 'grading', packet + '.grade')
        if os.path.isfile(os.path.join(g, 'grade.json')):
            shutil.copytree(g, os.path.join(dst, 'grade'))
            graded = True
        elif os.path.isdir(g):      # the grader answered but no attempt was accepted: kept under a name that says so
            shutil.copytree(g, os.path.join(dst, 'grade_not_accepted'))
    n_files = sum(len(fs) for _, _, fs in os.walk(os.path.join(dst, 'submission')))
    print('%s: %d submission files, %s messages, grade %s%s' % (
        name, n_files, n_msg, 'copied' if graded else 'none', '; not copied: %d' % len(notes) if notes else ''))
