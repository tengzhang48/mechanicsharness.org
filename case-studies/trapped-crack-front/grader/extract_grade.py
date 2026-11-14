"""Take the grade out of one grader run and decide whether it is accepted.

usage: python3 -I extract_grade.py OUT_JSONL PACKET_DIR GRADE_DIR

OUT_JSONL is the event stream of `opencode run --format json`. The grade is the last JSON value that follows a
```json fence in the grader's last text message. Files left in GRADE_DIR by an earlier attempt are removed first; then
  grade.md              the grader's whole last message
  run.json              tool calls, tokens and the result of the checks below
  grade.json            the grade, ONLY if it is accepted
  grade.rejected.json   the answer, if one was read but not accepted
A grade is accepted (exit 0) only if
  - its form is right: the allowed words, evidence given, the summary of question 2 follows from the listed checks;
  - every tool call that completed was read, glob or grep inside PACKET_DIR (a refused call is recorded, not counted);
  - the packet files are unchanged against PACKET_DIR.sha256.
Otherwise exit 1.

Another benchmark uses this tool with its own form: the environment variable GRADE_FORM names a JSON file with
  "allowed"         the allowed answers of the three questions (as ALLOWED below);
  "table"           the name of the table inside theory_source, and "table_keys" its rows (default: coefficients, COEF);
  "derive_summary"  true: the summary word of question 2 is derived here from the listed checks and written into the
                    grade; the grader's own word is kept as claimed_checks.answer_as_written and a difference is
                    recorded in run.json, not counted against the grade. Default false: a difference rejects the grade.
"""
import hashlib, json, os, sys

ALLOWED = {'theory_source': {'derivation_shown', 'own_computation', 'stated_with_source', 'stated_without_source', 'cannot_tell'},
           'claimed_checks': {'all_supported', 'some_unsupported', 'contradicted', 'cannot_tell', 'no_checks_claimed'},
           'error_location': {'none', 'theory', 'application', 'format', 'incomplete', 'cannot_tell'}}
STATUS = ('supported', 'unsupported', 'contradicted', 'cannot_tell')
COEF = ['p1', 'p2', 'p3', 'b2', 'b3', 'm3']
READING_TOOLS = ('read', 'glob', 'grep')
TABLE, DERIVE_SUMMARY = 'coefficients', False
if os.environ.get('GRADE_FORM'):
    _form = json.load(open(os.environ['GRADE_FORM']))
    ALLOWED = {k: set(v) for k, v in _form['allowed'].items()}
    TABLE, COEF, DERIVE_SUMMARY = _form['table'], list(_form['table_keys']), bool(_form.get('derive_summary'))
    assert set(ALLOWED) == {'theory_source', 'claimed_checks', 'error_location'}
OUTPUTS = ('grade.md', 'run.json', 'grade.json', 'grade.rejected.json')


def summary_of(statuses):
    """The summary word of question 2 that follows from the statuses of the listed checks."""
    if not statuses:
        return 'no_checks_claimed'
    for word, status in (('contradicted', 'contradicted'), ('some_unsupported', 'unsupported'), ('cannot_tell', 'cannot_tell')):
        if status in statuses:
            return word
    return 'all_supported'


def problems(g):
    out = []
    if not isinstance(g, dict):
        return ['the grade is not a JSON object']
    for key, allowed in ALLOWED.items():
        v = g.get(key)
        if key == 'claimed_checks' and DERIVE_SUMMARY and isinstance(v, dict):
            continue                      # the summary word is derived from the listed checks, whatever was written
        if not isinstance(v, dict) or v.get('answer') not in allowed:
            out.append('%s.answer is %r' % (key, v.get('answer') if isinstance(v, dict) else v))
    ts = g.get('theory_source') if isinstance(g.get('theory_source'), dict) else {}
    co = ts.get(TABLE) if isinstance(ts.get(TABLE), dict) else {}
    for c in COEF:
        if not isinstance(co.get(c), dict) or co[c].get('equals_key') not in ('yes', 'no', 'not_stated'):
            out.append('%s %s: equals_key missing or not yes/no/not_stated' % (TABLE, c))
    if not isinstance(ts.get('evidence'), list) or (not ts.get('evidence') and ts.get('answer') != 'cannot_tell'):
        out.append('theory_source.evidence is empty')
    cc = g.get('claimed_checks') if isinstance(g.get('claimed_checks'), dict) else {}
    checks = cc.get('checks')
    if not isinstance(checks, list):
        out.append('claimed_checks.checks is not a list')
    else:
        for i, c in enumerate(checks):
            if not isinstance(c, dict) or c.get('status') not in STATUS:
                out.append('claimed_checks.checks[%d].status missing or not allowed' % i)
        want = summary_of([c.get('status') for c in checks if isinstance(c, dict)])
        if cc.get('answer') in ALLOWED['claimed_checks'] and cc.get('answer') != want and not DERIVE_SUMMARY:
            out.append('claimed_checks.answer %r does not follow from the listed statuses (%s)' % (cc.get('answer'), want))
    el = g.get('error_location') if isinstance(g.get('error_location'), dict) else {}
    if el.get('answer') in ('theory', 'application', 'format', 'incomplete') and not str(el.get('first_wrong_step') or '').strip():
        out.append('error_location.first_wrong_step is empty')
    return out


def last_json(text):
    """(value,) for the last JSON value that follows a ```json fence, or None. The value is read by the JSON parser, so
    a fence inside one of its strings (a quoted code block) does not end it."""
    dec, found, pos = json.JSONDecoder(), None, 0
    while True:
        i = text.find('```json', pos)
        if i < 0:
            return found
        j = i + len('```json')
        while j < len(text) and text[j] in ' \t\r\n':
            j += 1
        try:
            value, pos = dec.raw_decode(text, j)
            found = (value,)
        except ValueError:
            pos = i + len('```json')


def leaves_packet(tool, inp, packet):
    """Why this tool input reaches outside the packet, or '' if it stays inside. Relative paths count from the packet,
    which is the folder the grader runs in."""
    def inside(path):
        p = os.path.realpath(os.path.join(packet, path))
        return p == packet or p.startswith(packet + os.sep)
    for key in ('path', 'filePath', 'file_path'):
        if isinstance(inp.get(key), str) and inp[key] and not inside(inp[key]):
            return '%s %s' % (key, inp[key])
    pattern = inp.get('pattern')
    if tool == 'glob' and isinstance(pattern, str) and (os.path.isabs(pattern) or '..' in pattern.replace('\\', '/').split('/')):
        fixed = pattern
        for ch in '*?[{':
            fixed = fixed.split(ch, 1)[0]
        base = inp.get('path') if isinstance(inp.get('path'), str) and inp.get('path') else packet
        if not inside(os.path.join(base, fixed)):
            return 'pattern %s' % pattern
    return ''


def main(stream, packet, grade_dir):
    packet = os.path.realpath(packet)
    os.makedirs(grade_dir, exist_ok=True)
    for name in OUTPUTS:                                    # nothing of an earlier attempt may pass for this one
        if os.path.lexists(os.path.join(grade_dir, name)):
            os.remove(os.path.join(grade_dir, name))
    texts, tokens, errors = [], [], []
    summary_note = ''
    tools, paths_read, refused, violations = [], set(), [], []
    for line in open(stream, encoding='utf-8', errors='replace'):
        try:
            j = json.loads(line)
        except ValueError:
            continue
        if not isinstance(j, dict):
            continue
        part = j.get('part') if isinstance(j.get('part'), dict) else {}
        if j.get('type') == 'text' and part.get('text'):
            texts.append(part['text'])
        elif j.get('type') == 'tool_use':
            st = part.get('state') if isinstance(part.get('state'), dict) else {}
            inp = st.get('input') if isinstance(st.get('input'), dict) else {}
            tool, shown = part.get('tool'), json.dumps(st.get('input'), ensure_ascii=False)[:300]
            tools.append(tool)
            if st.get('status') != 'completed':
                refused.append({'tool': tool, 'status': st.get('status'), 'input': shown,
                                'error': str(st.get('error') or st.get('output') or '')[:200]})
            elif tool not in READING_TOOLS:
                violations.append('the grader used the tool %r (input %s)' % (tool, shown[:120]))
            else:
                why = leaves_packet(tool, inp, packet)
                if why:
                    violations.append('%s outside the packet: %s' % (tool, why))
                elif tool == 'read' and isinstance(inp.get('path'), str):
                    paths_read.add(inp['path'])
        elif j.get('type') == 'step_finish' and isinstance(part.get('tokens'), dict):
            tokens.append(part['tokens'])
        elif j.get('type') == 'error':
            errors.append(json.dumps(j)[:500])
    last = texts[-1] if texts else ''
    with open(os.path.join(grade_dir, 'grade.md'), 'w', encoding='utf-8') as fh:
        fh.write(last + '\n')
    grade, bad = None, []
    found = last_json(last)
    if found is None:
        bad.append('the json block does not parse' if '```json' in last else 'no fenced json block in the last message')
    else:
        grade = found[0]
        bad += problems(grade)
        cc = grade.get('claimed_checks') if isinstance(grade, dict) else None
        if DERIVE_SUMMARY and isinstance(cc, dict) and isinstance(cc.get('checks'), list):
            want = summary_of([c.get('status') for c in cc['checks'] if isinstance(c, dict)])
            if cc.get('answer') != want:
                summary_note = 'summary of question 2 written by the grader: %r; derived from its listed checks: %r' % (cc.get('answer'), want)
            cc['answer_as_written'], cc['answer'] = cc.get('answer'), want
    bad += violations
    changed, listed = [], {}
    for row in open(packet + '.sha256'):
        h, rel = row.rstrip('\n').split('  ', 1)
        listed[rel] = h
    now = {os.path.relpath(os.path.join(d, f), packet) for d, _, fs in os.walk(packet) for f in fs}
    for rel in sorted(now | set(listed)):
        p = os.path.join(packet, rel)
        if rel not in listed:
            changed.append('new file ' + rel)
        elif not os.path.isfile(p) or os.path.islink(p):
            changed.append('missing ' + rel)
        else:
            with open(p, 'rb') as fh:
                if hashlib.sha256(fh.read()).hexdigest() != listed[rel]:
                    changed.append('changed ' + rel)
    if changed:
        bad.append('the packet changed during grading: ' + '; '.join(changed))
    run = {'accepted': not bad, 'form_problems': bad, 'tool_calls': len(tools), 'tools_used': sorted({str(t) for t in tools}),
           'paths_read': sorted(paths_read), 'refused_tool_calls': refused, 'packet_unchanged': not changed,
           'steps': len(tokens), 'reasoning_tokens': sum(t.get('reasoning', 0) for t in tokens),
           'input_tokens': sum(t.get('input', 0) for t in tokens), 'output_tokens': sum(t.get('output', 0) for t in tokens),
           'stream_errors': errors}
    if summary_note:
        run['summary_note'] = summary_note
    with open(os.path.join(grade_dir, 'run.json'), 'w') as fh:
        json.dump(run, fh, indent=1)
    if grade is not None:
        with open(os.path.join(grade_dir, 'grade.rejected.json' if bad else 'grade.json'), 'w') as fh:
            json.dump(grade, fh, indent=1, ensure_ascii=False)
    if bad:
        print('GRADE NOT ACCEPTED:', '; '.join(bad))
        return 1
    print('grade: theory_source %s | claimed_checks %s | error_location %s | %d tool calls (%d refused), %d reasoning tokens' %
          (grade['theory_source']['answer'], grade['claimed_checks']['answer'], grade['error_location']['answer'],
           run['tool_calls'], len(refused), run['reasoning_tokens']))
    return 0


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], os.path.abspath(sys.argv[2]), sys.argv[3]))
