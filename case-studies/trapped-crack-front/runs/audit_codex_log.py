"""Audit of one Codex run after it ended, from its text log: the shell commands it ran and any sign that it reached the
network other than through the model. Organizer only. The counterpart of audit_log.py for the Codex CLI.

usage: python3 -I audit_codex_log.py RUN_DIR [messages]
  without the second word: prints the audit (the report kept as AUDIT.md);
  with it: prints the agent's visible messages and commands in order (kept as messages.md).

The log (RUN_DIR/codex.log) is a sequence of blocks that begin with a line "user", "codex" (a message of the agent),
"exec" (a shell command, then its result) or "thinking" (a summary of reasoning, if the run printed any). Reported:
the number of lines, the tokens the log states, every command (its first line), every command or message line that
mentions the network (the patterns of audit_log.py) or a path outside /work. A hit is a line for a person to read,
not a verdict. The output of commands is not searched: it would repeat the task text and the agent's own files.
"""
import os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_log import NETWORK, OUTSIDE   # noqa: E402

HEADS = ('user', 'codex', 'exec', 'thinking', 'tokens used')
RESULT = re.compile(r'^ (succeeded|exited \d+|failed|timed out).* in [\d.]+m?s:?')


def blocks(path):
    """[(kind, [lines])] in order. Inside an exec block, lines after the result line are the command's output."""
    out, kind, cur = [], None, []
    for raw in open(path, encoding='utf-8', errors='replace'):
        line = raw.rstrip('\n')
        if line in HEADS:
            if kind:
                out.append((kind, cur))
            kind, cur = line, []
        elif kind:
            cur.append(line)
    if kind:
        out.append((kind, cur))
    return out


def command_of(lines):
    """(the command text, its result line) of an exec block."""
    cmd = []
    for ln in lines:
        if RESULT.match(ln):
            return '\n'.join(cmd), ln.strip()
        cmd.append(ln)
    return '\n'.join(cmd), 'no result line'


def short(x, n=200):
    x = ' '.join(x.split())
    return x if len(x) <= n else x[:n] + ' ...'


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    run = sys.argv[1]
    name = os.path.basename(os.path.normpath(run))
    log = os.path.join(run, 'codex.log')
    bl = blocks(log)
    n_lines = sum(1 for _ in open(log, encoding='utf-8', errors='replace'))
    if len(sys.argv) == 3:
        print('# Visible messages and commands of %s, in order\n' % name)
        print('Cut from the log of the Codex CLI (codex.log, kept outside git). Commands are shortened; their output is left out.\n')
        k = 0
        for kind, lines in bl:
            if kind == 'codex':
                k += 1
                print('\n## Message %d\n\n%s\n' % (k, '\n'.join(lines).strip()))
            elif kind == 'thinking':
                print('\n*Reasoning summary:* %s\n' % short('\n'.join(lines), 600))
            elif kind == 'exec':
                cmd, res = command_of(lines)
                print('- `%s` (%s)' % (short(cmd, 300).replace('`', "'"), res))
        sys.exit(0)
    cmds = [command_of(lines) for kind, lines in bl if kind == 'exec']
    said = ['\n'.join(lines) for kind, lines in bl if kind in ('codex', 'thinking')]
    tokens = [short(' '.join(lines), 40) for kind, lines in bl if kind == 'tokens used']
    net, outside = [], []
    for what, text in [('command', c) for c, _ in cmds] + [('message', m) for m in said]:
        for ln in text.split('\n'):
            hits = sorted(set(m.group(0).lower() for m in NETWORK.finditer(ln)))
            if hits:
                net.append('%s [%s]: %s' % (what, ', '.join(hits), short(ln, 220)))
            if what == 'command':
                outs = sorted(set(m.group(0) for m in OUTSIDE.finditer(ln)))
                outs = [o for o in outs if not o.startswith(('/usr/bin/bash', '/usr/bin/env', '/usr/bin/python'))]
                if outs:
                    outside.append('%s: %s' % (', '.join(outs), short(ln, 200)))
    print('# Audit of %s\n' % name)
    print('Codex CLI log (codex.log, %d lines). Tokens used, as the log states: %s.\n' % (n_lines, tokens[-1].split()[0] if tokens else 'not stated'))
    print('Commands: %d; of these %d did not succeed (a failed command is usually a program error, listed below).\n' % (
        len(cmds), sum(1 for _, r in cmds if not r.startswith('succeeded'))))
    print('## Shell commands (%d)\n' % len(cmds))
    for c, r in cmds:
        print('- %s  [%s]' % (short(c.split('\n')[0], 220), r))
    print('\n## Lines of commands or messages that mention the network, a web tool or a connector (%d)\n' % len(net))
    for x in net:
        print('- ' + x)
    print('\n## Commands that name a path outside /work (%d)\n' % len(outside))
    for x in outside:
        print('- ' + x)
    print('\nSettings of the run: web_search = "disabled", [features] apps = false (see launch_codex.sh).')
