"""Assemble the folder one grader run reads: the grading instructions, the task as given, the key notes, the scorer's
report and a copy of the solver's folder. Organizer only.

usage: python3 -I make_packet.py SUBMISSION_DIR SCORE_DIR PACKET_DIR [TASK_FILE]      (PACKET_DIR must not exist)

- TASK_FILE is the task text as the solver received it. Default: ../task/PUBLIC_TASK.md. For a run of a task variant
  (for example the reminder variant) give the PUBLIC_TASK.md of that run's workspace, so that TASK.md in the packet is
  word for word what the solver saw.

- Give PACKET_DIR a neutral name (a run number, not a model name): the grader should not know who solved the problem.
- Only regular files that are UTF-8 text of at most 300 kB are copied from the submission. Everything else is listed in
  MANIFEST.txt with its size and the reason it was left out: larger files, files that are not UTF-8 text, symbolic
  links (never followed, whatever they point to) and the folders of SKIP_DIRS.
- A folder without results.json is packed too, so that a run that did not finish can still be read; the manifest says
  that the file is absent. The scorer's report (SCORE_DIR/score.md) is required.
- Writes PACKET_DIR.sha256 beside the packet, to show afterwards that the grader changed nothing.
- Another benchmark uses this tool with its own instructions: the environment variable GRADER_DIR names the folder
  that holds its GRADING.md and KEY_NOTES.md (default: this folder); the default TASK_FILE is then
  GRADER_DIR/../task/PUBLIC_TASK.md.
- PACKET_BUNDLE=1 (environment, optional) also writes SUBMISSION_ALL.txt into the packet: every copied file of the
  submission in one file, each line with its line number, so that the grader needs a few tool calls instead of one
  per file (opencode's "Duplicate tool input start" error grows with the number of calls). A copy of the task text
  in the submission that is identical to TASK.md is named there, not repeated.
"""
import hashlib, os, shutil, stat, sys

HERE = os.path.abspath(os.environ.get('GRADER_DIR') or os.path.dirname(os.path.abspath(__file__)))
LIMIT = 300_000
SKIP_DIRS = {'.git', '__pycache__', 'node_modules', '.venv', 'venv', '.cache', '.config', '.local'}


def is_text(path):
    """UTF-8 text without NUL bytes. Called only for files of at most LIMIT bytes, so the whole file is read."""
    with open(path, 'rb') as fh:
        data = fh.read()
    if b'\0' in data:
        return False
    try:
        data.decode('utf-8')
        return True
    except UnicodeDecodeError:
        return False


def main(sub, score_dir, packet, task_file=None):
    task_file = task_file or os.path.join(HERE, '..', 'task', 'PUBLIC_TASK.md')
    if not os.path.isfile(task_file):
        sys.exit('missing task text %s' % task_file)
    if os.path.lexists(packet):
        sys.exit('%s exists; refusing to overwrite' % packet)
    if not os.path.isdir(sub):
        sys.exit('missing folder %s' % sub)
    if not os.path.isfile(os.path.join(score_dir, 'score.md')):
        sys.exit('missing %s' % os.path.join(score_dir, 'score.md'))
    os.makedirs(os.path.join(packet, 'submission'))
    shutil.copyfile(os.path.join(HERE, 'GRADING.md'), os.path.join(packet, 'GRADING.md'))
    shutil.copyfile(os.path.join(HERE, 'KEY_NOTES.md'), os.path.join(packet, 'KEY_NOTES.md'))
    shutil.copyfile(task_file, os.path.join(packet, 'TASK.md'))
    shutil.copyfile(os.path.join(score_dir, 'score.md'), os.path.join(packet, 'score.md'))
    lines = []

    def left_out(size, why, rel):
        lines.append('%10d  LEFT OUT (%s)  %s' % (size, why, rel))

    for d, dirs, files in os.walk(sub):                     # os.walk does not go into folders that are links
        keep = []
        for name in sorted(dirs):
            p = os.path.join(d, name)
            if os.path.islink(p):
                left_out(os.lstat(p).st_size, 'symbolic link', os.path.relpath(p, sub))
            elif name in SKIP_DIRS:
                left_out(0, 'folder not copied', os.path.relpath(p, sub))
            else:
                keep.append(name)
        dirs[:] = keep
        for name in sorted(files):
            src = os.path.join(d, name)
            rel = os.path.relpath(src, sub)
            st = os.lstat(src)                              # the entry itself: a link is not followed
            if stat.S_ISLNK(st.st_mode):
                left_out(st.st_size, 'symbolic link', rel)
            elif not stat.S_ISREG(st.st_mode):
                left_out(0, 'not a regular file', rel)
            elif st.st_size > LIMIT:
                left_out(st.st_size, 'too large', rel)
            elif not is_text(src):
                left_out(st.st_size, 'not a text file', rel)
            else:
                dst = os.path.join(packet, 'submission', rel)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copyfile(src, dst)
                lines.append('%10d  copied  %s' % (st.st_size, rel))
    if os.environ.get('PACKET_BUNDLE') == '1':
        task_bytes = open(task_file, 'rb').read()
        with open(os.path.join(packet, 'SUBMISSION_ALL.txt'), 'w', encoding='utf-8') as fh:
            fh.write('Every copied file of submission/, one after the other. Each line begins with its line number in that file.\n')
            for row in lines:
                if '  copied  ' not in row:
                    continue
                rel = row.split('  copied  ', 1)[1]
                data = open(os.path.join(packet, 'submission', rel), 'rb').read()
                if data == task_bytes:
                    fh.write('\n===== submission/%s: identical to TASK.md, not repeated here =====\n' % rel)
                    continue
                text = data.decode('utf-8').split('\n')
                if text and text[-1] == '':
                    text.pop()
                fh.write('\n===== submission/%s (%d lines) =====\n' % (rel, len(text)))
                fh.write(''.join('%5d  %s\n' % (i + 1, ln) for i, ln in enumerate(text)))
    has_results = os.path.isfile(os.path.join(packet, 'submission', 'results.json'))
    with open(os.path.join(packet, 'MANIFEST.txt'), 'w') as fh:
        fh.write('Files of the solver\'s folder (bytes, status, path under submission/)\n')
        fh.write('results.json is present.\n' if has_results else
                 'results.json is ABSENT: the copy of the solver\'s folder has no results.json at its top level.\n')
        fh.write('\n'.join(lines) + '\n')
    sums = []
    for d, _, files in os.walk(packet):
        for name in sorted(files):
            p = os.path.join(d, name)
            with open(p, 'rb') as fh:
                sums.append('%s  %s' % (hashlib.sha256(fh.read()).hexdigest(), os.path.relpath(p, packet)))
    with open(packet.rstrip('/') + '.sha256', 'w') as fh:
        fh.write('\n'.join(sorted(sums, key=lambda s: s.split('  ', 1)[1])) + '\n')
    print('packet', packet, 'with', len(lines), 'entries listed,', sum('  copied  ' in x for x in lines), 'copied;',
          'results.json present' if has_results else 'results.json ABSENT')


if __name__ == '__main__':
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    main(*[os.path.abspath(a) for a in sys.argv[1:]])
