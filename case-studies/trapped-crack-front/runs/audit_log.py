"""Audit of one opencode run after it ended: what the agent did with its tools, and any sign that it reached the
network other than through the model. Organizer only.

usage: python3 -I audit_log.py RUN_DIR        (reads RUN_DIR/opencode.jsonl, prints the report)

Reported: the tool calls by tool and outcome; every call that did not complete (a refused web request is one); every
call whose input mentions the network (a URL, a network library or program, opencode's web tools), with the outcome;
every path outside /work that a file tool was given and every shell command that names a folder outside /work; the
tokens and the time span of the stream. A hit is a line for a person to read, not a verdict: file contents written
by the agent are searched too, so a word such as "socket" in a comment shows up here.
"""
import collections, json, os, re, sys

NETWORK = re.compile(r'https?://|ftp://|\burllib\b|\brequests\b|\bhttpx\b|\baiohttp\b|\bhttp\.client\b|\bsocket\b|'
                     r'\bcurl\b|\bwget\b|\bpip3?\s+install\b|\bapt(?:-get)?\s+install\b|\bgit\s+clone\b|\bfetch\s*\(|'
                     r'\bwebfetch\b|\bwebsearch\b|\bssh\b|\bscp\b|\bnc\s+-|\barxiv\b|\bdoi\.org\b|\bscholar\b', re.I)
OUTSIDE = re.compile(r'(?<![\w./-])/(?:ochome|home|root|etc|proc|sys|opt|media|mnt|var|usr/local|srv)\b[^\s"\']*')
PATH_KEYS = ('path', 'filePath', 'file_path', 'cwd', 'workdir')


def main(run):
    tools, not_done, net, outside = collections.Counter(), [], [], []
    tokens = collections.Counter()
    first = last = None
    texts, errors, steps = [], [], 0
    for line in open(os.path.join(run, 'opencode.jsonl'), encoding='utf-8', errors='replace'):
        try:
            j = json.loads(line)
        except ValueError:
            continue
        if not isinstance(j, dict):
            continue
        t = j.get('timestamp')
        if isinstance(t, (int, float)):
            first, last = (t if first is None else first), t
        part = j.get('part') if isinstance(j.get('part'), dict) else {}
        kind = j.get('type')
        if kind == 'text' and part.get('text'):
            texts.append(part['text'])
        elif kind == 'error':
            errors.append(json.dumps(j.get('error'))[:300])
        elif kind == 'step_finish':
            steps += 1
            for k, v in (part.get('tokens') or {}).items():
                if isinstance(v, (int, float)):
                    tokens[k] += v
                elif isinstance(v, dict):
                    for kk, vv in v.items():
                        if isinstance(vv, (int, float)):
                            tokens[k + '.' + kk] += vv
        elif kind == 'tool_use':
            st = part.get('state') if isinstance(part.get('state'), dict) else {}
            inp = st.get('input') if isinstance(st.get('input'), dict) else {}
            tool, status = str(part.get('tool')), str(st.get('status'))
            tools[(tool, status)] += 1
            shown = json.dumps(inp, ensure_ascii=False)
            result = str(st.get('error') or st.get('output') or '')
            if status != 'completed':
                not_done.append((tool, status, shown[:200], result[:160]))
            hits = sorted(set(m.group(0).lower() for m in NETWORK.finditer(shown)))
            if hits or tool in ('webfetch', 'websearch'):
                net.append((tool, status, hits, shown[:300], result[:200]))
            for k in PATH_KEYS:
                v = inp.get(k)
                if isinstance(v, str) and v.startswith('/') and not (v == '/work' or v.startswith('/work/')):
                    outside.append((tool, status, '%s = %s' % (k, v), result[:120]))
            if tool not in ('read', 'write', 'edit', 'glob', 'grep'):
                for m in OUTSIDE.finditer(' '.join(str(v) for k, v in inp.items() if k not in ('content',))):
                    outside.append((tool, status, m.group(0), shown[:200]))
    print('# Audit of %s' % os.path.basename(os.path.abspath(run)))
    if first is not None:
        print('\nStream: %d steps, %.1f minutes between the first and the last event.' % (steps, (last - first) / 60000.0))
    print('Tokens: ' + ', '.join('%s %d' % kv for kv in sorted(tokens.items())))
    print('Errors in the stream: %d%s' % (len(errors), ''.join('\n  ' + e for e in errors[:5])))
    print('\n## Tool calls\n')
    for (tool, status), n in sorted(tools.items()):
        print('- %s, %s: %d' % (tool, status, n))
    print('\n## Calls that did not complete (%d)\n' % len(not_done))
    for row in not_done:
        print('- %s (%s): %s -> %s' % row)
    print('\n## Calls whose input mentions the network (%d)\n' % len(net))
    for tool, status, hits, shown, result in net:
        print('- %s (%s), matched %s: %s -> %s' % (tool, status, hits, shown, result.replace('\n', ' ')))
    print('\n## Paths outside /work (%d)\n' % len(outside))
    for row in outside:
        print('- %s (%s): %s -> %s' % tuple(str(x).replace('\n', ' ') for x in row))
    print('\n## The agent\'s last message\n')
    print(texts[-1] if texts else '(none)')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
