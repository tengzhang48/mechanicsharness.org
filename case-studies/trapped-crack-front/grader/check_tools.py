"""Tests of the grader's two tools on made-up inputs: make_packet.py (what is copied for the grader) and
extract_grade.py (what is accepted as a grade). No model is called. Organizer only.

usage: python3 -I check_tools.py          exit 1 if a test fails

Most of these cases come from the review of draft 0.2 by Gemini (13 Nov 2026): files in a solver's folder that are
not plain text, links, a missing results.json, a long tool input, a rejected grade left on disk, a fenced block inside
a quotation, and tool calls that leave the packet.
"""
import copy, gzip, hashlib, json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
results = []


def check(name, ok, detail=''):
    results.append(bool(ok))
    print('%s  %s%s' % ('ok  ' if ok else 'FAIL', name, '   [%s]' % detail if detail and not ok else ''))


def run(script, *args):
    r = subprocess.run([sys.executable, '-I', os.path.join(HERE, script)] + list(args), capture_output=True, text=True)
    return r.returncode, (r.stdout + r.stderr).strip()


def write(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb' if isinstance(data, bytes) else 'w') as fh:
        fh.write(data)


def manifest_line(packet, rel):
    for line in open(os.path.join(packet, 'MANIFEST.txt')).read().splitlines():
        if line.endswith('  ' + rel):
            return line
    return ''


def left_out(packet, rel):
    return 'LEFT OUT' in manifest_line(packet, rel) and not os.path.lexists(os.path.join(packet, 'submission', rel))


def copied(packet, rel, src):
    dst = os.path.join(packet, 'submission', rel)
    return ('copied' in manifest_line(packet, rel) and os.path.isfile(dst) and not os.path.islink(dst)
            and open(dst, 'rb').read() == open(src, 'rb').read())


# ---------------------------------------------------------------------------------------------- make_packet.py
def packet_tests(tmp):
    sub, score, outside = (os.path.join(tmp, x) for x in ('sub', 'score', 'outside'))
    write(os.path.join(score, 'score.md'), '# score\n')
    write(os.path.join(outside, 'secret.txt'), 'SECRET-OUTSIDE\n')
    write(os.path.join(outside, 'dir', 'inner.txt'), 'SECRET-INNER\n')
    write(os.path.join(sub, 'results.json'), '{}\n')
    write(os.path.join(sub, 'method.md'), '# Method\n')
    write(os.path.join(sub, 'deep', 'er', 'code.py'), 'print(1)\n')
    write(os.path.join(sub, 'big.txt'), 'a' * 300_001)
    write(os.path.join(sub, 'bytes_no_nul.bin'), bytes([0xE9, 0x80, 0xFF, 0xC3]) * 1500)      # 6000 bytes, no NUL, not UTF-8
    write(os.path.join(sub, 'late_nul.dat'), b'a' * 5000 + b'\0' + b'a' * 10)                 # NUL after the first 4 kB
    write(os.path.join(sub, 'utf8_at_4096.txt'), b'a' * 4095 + 'é'.encode() + b'b' * 100)     # a 2-byte character across byte 4096
    write(os.path.join(sub, 'latin1.txt'), b'caf\xe9\n')
    write(os.path.join(sub, 'data.gz'), gzip.compress(b'hello ' * 2000, mtime=0))
    write(os.path.join(sub, '.git', 'config'), '[core]\n')
    write(os.path.join(sub, '__pycache__', 'x.pyc'), b'\0\1\2')
    os.symlink('nowhere', os.path.join(sub, 'broken_link'))
    os.symlink(os.path.join(outside, 'secret.txt'), os.path.join(sub, 'link_to_file.txt'))
    os.symlink(os.path.join(outside, 'dir'), os.path.join(sub, 'link_to_dir'))

    packet = os.path.join(tmp, 'P1')
    code, out = run('make_packet.py', sub, score, packet)
    check('packet: a folder with a broken link, links and binary files is packed without a crash', code == 0, out[-300:])
    if code != 0:
        for name in ('the rest of the packet tests',):
            check('packet: %s (not run: no packet)' % name, False)
        return None
    for rel in ('results.json', 'method.md', 'deep/er/code.py', 'utf8_at_4096.txt'):
        check('packet: %s is copied unchanged' % rel, copied(packet, rel, os.path.join(sub, rel)), manifest_line(packet, rel))
    for rel, why in (('big.txt', 'larger than 300 kB'), ('bytes_no_nul.bin', '6000 bytes that are not UTF-8, no NUL'),
                     ('late_nul.dat', 'a NUL byte after the first 4 kB'), ('latin1.txt', 'not UTF-8'),
                     ('data.gz', 'a gzip stream'), ('broken_link', 'a link to nothing'),
                     ('link_to_file.txt', 'a link to a file outside the folder'),
                     ('link_to_dir', 'a link to a folder outside')):
        check('packet: %s (%s) is listed as left out and not copied' % (rel, why), left_out(packet, rel), manifest_line(packet, rel))
    for rel in ('.git', '__pycache__'):
        check('packet: the folder %s is listed as left out and not copied' % rel, left_out(packet, rel), manifest_line(packet, rel))
    everything = b''.join(open(os.path.join(d, f), 'rb').read() for d, _, fs in os.walk(packet) for f in fs)
    check('packet: nothing from outside the solver\'s folder is in the packet', b'SECRET' not in everything)
    listed = dict(reversed(row.split('  ', 1)) for row in open(packet + '.sha256').read().splitlines())
    now = {os.path.relpath(os.path.join(d, f), packet) for d, _, fs in os.walk(packet) for f in fs}
    check('packet: the hash list names every file of the packet, with the right hash',
          set(listed) == now and all(hashlib.sha256(open(os.path.join(packet, r), 'rb').read()).hexdigest() == h
                                     for r, h in listed.items()))
    code, out = run('make_packet.py', sub, score, packet)
    check('packet: an existing packet folder is not overwritten', code != 0)

    sub2 = os.path.join(tmp, 'sub2')                        # a run that wrote its method and code but no results.json
    write(os.path.join(sub2, 'method.md'), '# Method\nnot finished\n')
    write(os.path.join(sub2, 'out', 'result.json'), '{}\n')
    p2 = os.path.join(tmp, 'P2')
    code, out = run('make_packet.py', sub2, score, p2)
    check('packet: a folder without results.json is still packed', code == 0 and os.path.isfile(os.path.join(p2, 'submission', 'method.md')), out[-200:])
    check('packet: its manifest says that results.json is absent',
          code == 0 and 'results.json' in open(os.path.join(p2, 'MANIFEST.txt')).read().split('\n')[1])
    variant = os.path.join(tmp, 'variant_task.md')            # a run of a task variant: TASK.md is the text that run received
    write(variant, 'the task with one more paragraph\n')
    p4 = os.path.join(tmp, 'P4')
    code, out = run('make_packet.py', sub2, score, p4, variant)
    check('packet: with a task file given, TASK.md is that file', code == 0 and open(os.path.join(p4, 'TASK.md')).read() == open(variant).read(), out[-200:])
    check('packet: without one, TASK.md is the task text of the benchmark',
          open(os.path.join(p2, 'TASK.md')).read() == open(os.path.join(HERE, '..', 'task', 'PUBLIC_TASK.md')).read())
    code, out = run('make_packet.py', sub2, score, os.path.join(tmp, 'P5'), os.path.join(tmp, 'no_such_task.md'))
    check('packet: no packet with a task file that does not exist', code != 0 and not os.path.exists(os.path.join(tmp, 'P5')))
    code, out = run('make_packet.py', sub2, os.path.join(tmp, 'no_score'), os.path.join(tmp, 'P3'))
    check('packet: no packet without the scorer\'s report', code != 0 and not os.path.exists(os.path.join(tmp, 'P3')))
    return packet


# -------------------------------------------------------------------------------------------- extract_grade.py
GOOD = {
    'theory_source': {'answer': 'stated_with_source',
                      'evidence': [{'file': 'submission/method.md', 'lines': '1-2', 'quote': 'Gao and Rice (1986)'}],
                      'coefficients': {c: {'used': 'x', 'equals_key': 'yes'} for c in ('p1', 'p2', 'p3', 'b2', 'b3', 'm3')}},
    'claimed_checks': {'answer': 'all_supported',
                       'checks': [{'claim': 'G = Gamma', 'claimed_at': 'submission/method.md:9', 'status': 'supported',
                                   'evidence': 'submission/solve.py:150'}]},
    'error_location': {'answer': 'none', 'first_wrong_step': '', 'where': '', 'explains_groups': []},
    'notes': ''}


def tool(name, inp, status='completed', output='text of the file'):
    return {'type': 'tool_use', 'part': {'tool': name, 'state': {'status': status, 'input': inp, 'output': output}}}


def answer(grade, before='I read every file.\n\n'):
    return {'type': 'text', 'part': {'text': before + '```json\n' + json.dumps(grade, indent=1) + '\n```\n'}}


STEP = {'type': 'step_finish', 'part': {'tokens': {'input': 100, 'output': 50, 'reasoning': 7}}}


def extract(tmp, packet, name, events, grade_dir=None):
    path = os.path.join(tmp, name + '.jsonl')
    with open(path, 'w') as fh:
        fh.write('\n'.join(json.dumps(e) for e in events) + '\n')
    gd = grade_dir or os.path.join(tmp, name + '.grade')
    code, out = run('extract_grade.py', path, packet, gd)
    g = os.path.join(gd, 'grade.json')
    r = os.path.join(gd, 'run.json')
    return code, out, (json.load(open(g)) if os.path.isfile(g) else None), (json.load(open(r)) if os.path.isfile(r) else None), gd


def changed(grade, path, value):
    g = copy.deepcopy(grade)
    node = g
    for k in path[:-1]:
        node = node[k]
    node[path[-1]] = value
    return g


def grade_tests(tmp, packet):
    reads = [tool('read', {'path': os.path.join(packet, f)}) for f in ('GRADING.md', 'submission/method.md')]
    reads.append(tool('glob', {'pattern': '**/*'}))
    reads.append(tool('grep', {'pattern': 'Gamma', 'path': os.path.join(packet, 'submission')}))

    code, out, g, r, _ = extract(tmp, packet, 'e_good', reads + [STEP, answer(GOOD)])
    check('grade: a well-formed grade is accepted and written as given', code == 0 and g == GOOD and r and r['packet_unchanged'], out[-300:])

    long_path = os.path.join(packet, 'submission', 'd' * 150, 'a_long_file_name_for_a_check.py')
    code, out, g, r, _ = extract(tmp, packet, 'e_long', reads + [tool('read', {'path': long_path, 'offset': 1, 'limit': 2000}), STEP, answer(GOOD)])
    check('grade: a read whose input is longer than 200 characters does not stop the extraction',
          code == 0 and g == GOOD and r is not None and long_path in r.get('paths_read', []), out[-300:])

    bad_word = changed(GOOD, ('claimed_checks', 'answer'), 'supported')
    code, out, g, r, gd = extract(tmp, packet, 'e_word', reads + [STEP, answer(bad_word)])
    check('grade: a summary word that is not allowed is rejected, and no grade.json is written',
          code != 0 and g is None and r is not None and r['form_problems'], 'exit %d, grade.json %s' % (code, 'present' if g else 'absent'))
    check('grade: the rejected answer is kept under another name', os.path.isfile(os.path.join(gd, 'grade.rejected.json')))

    gd = os.path.join(tmp, 'e_stale.grade')                 # attempt 1 accepted, attempt 2 ends with an error and no answer
    extract(tmp, packet, 'e_stale_1', reads + [STEP, answer(GOOD)], gd)
    code, out, g, r, _ = extract(tmp, packet, 'e_stale_2', reads[:1] + [{'type': 'error', 'error': {'message': 'Duplicate tool input start'}}], gd)
    check('grade: a run without an answer leaves no grade.json from an earlier attempt', code != 0 and g is None, 'exit %d' % code)

    quoting = changed(GOOD, ('theory_source', 'evidence'), [{'file': 'submission/method.md', 'lines': '3-5',
                                                              'quote': '```python\nres = 1e-16\n```'}])
    code, out, g, r, _ = extract(tmp, packet, 'e_fence', reads + [STEP, answer(quoting)])
    check('grade: a fenced block quoted inside the grade does not break it', code == 0 and g == quoting, out[-300:])

    code, out, g, r, _ = extract(tmp, packet, 'e_twice', reads + [STEP, answer(bad_word, before='A first try:\n\n```json\n{"draft": 1}\n```\n\nFinal:\n\n')])
    check('grade: the last json block is the grade', code != 0 and r is not None and 'supported' in json.dumps(r['form_problems']), out[-300:])

    outside = os.path.join(os.path.dirname(packet), 'MAPPING.json')
    code, out, g, r, _ = extract(tmp, packet, 'e_outside', reads + [tool('read', {'path': outside}), STEP, answer(GOOD)])
    check('grade: a completed read outside the packet rejects the grade', code != 0 and g is None, out[-300:])
    code, out, g, r, _ = extract(tmp, packet, 'e_glob_up', reads + [tool('glob', {'pattern': '../*'}), STEP, answer(GOOD)])
    check('grade: a completed file search above the packet rejects the grade', code != 0 and g is None, out[-300:])
    code, out, g, r, _ = extract(tmp, packet, 'e_execute', reads + [tool('execute', {'code': 'return 17*23'}, output='391'), STEP, answer(GOOD)])
    check('grade: a completed call of a tool other than read, glob, grep rejects the grade', code != 0 and g is None, out[-300:])
    code, out, g, r, _ = extract(tmp, packet, 'e_refused', reads + [tool('read', {'path': outside}, status='error', output=None), STEP, answer(GOOD)])
    check('grade: a refused read outside the packet is recorded, and the grade stands',
          code == 0 and g == GOOD and r and len(r.get('refused_tool_calls', [])) == 1, out[-300:])

    one_unknown = changed(changed(GOOD, ('claimed_checks', 'checks'), GOOD['claimed_checks']['checks'] + [
        {'claim': 'Nystrom', 'claimed_at': 'submission/method.md:12', 'status': 'cannot_tell', 'evidence': 'output file left out'}]),
        ('claimed_checks', 'answer'), 'cannot_tell')
    code, out, g, r, _ = extract(tmp, packet, 'e_ct', reads + [STEP, answer(one_unknown)])
    check('grade: a check marked cannot_tell is accepted when the summary says cannot_tell', code == 0 and g == one_unknown, out[-300:])
    code, out, g, r, _ = extract(tmp, packet, 'e_ct2', reads + [STEP, answer(changed(one_unknown, ('claimed_checks', 'answer'), 'all_supported'))])
    check('grade: a summary that does not follow from the listed checks is rejected', code != 0 and g is None, out[-300:])
    no_source = changed(changed(GOOD, ('theory_source', 'answer'), 'cannot_tell'), ('theory_source', 'evidence'), [])
    code, out, g, r, _ = extract(tmp, packet, 'e_ct3', reads + [STEP, answer(no_source)])
    check('grade: theory_source cannot_tell needs no evidence', code == 0 and g == no_source, out[-300:])
    code, out, g, r, _ = extract(tmp, packet, 'e_noev', reads + [STEP, answer(changed(GOOD, ('theory_source', 'evidence'), []))])
    check('grade: any other theory_source answer without evidence is rejected', code != 0 and g is None, out[-300:])

    unfinished = changed(GOOD, ('error_location',), {'answer': 'incomplete', 'first_wrong_step': 'Part B was not computed.',
                                                     'where': 'submission/method.md:20', 'explains_groups': ['Part B']})
    code, out, g, r, _ = extract(tmp, packet, 'e_inc', reads + [STEP, answer(unfinished)])
    check('grade: error_location incomplete is accepted', code == 0 and g == unfinished, out[-300:])
    code, out, g, r, _ = extract(tmp, packet, 'e_inc2', reads + [STEP, answer(changed(unfinished, ('error_location', 'first_wrong_step'), ''))])
    check('grade: an error without a named step is rejected', code != 0 and g is None, out[-300:])
    code, out, g, r, _ = extract(tmp, packet, 'e_none', reads + [STEP, {'type': 'text', 'part': {'text': 'No grade.'}}])
    check('grade: an answer without a json block is rejected', code != 0 and g is None, out[-300:])

    write(os.path.join(packet, 'submission', 'method.md'), '# Method (changed)\n')          # last: these change the packet
    code, out, g, r, _ = extract(tmp, packet, 'e_changed', reads + [STEP, answer(GOOD)])
    check('grade: a packet file changed during grading rejects the grade', code != 0 and g is None and r and not r['packet_unchanged'], out[-300:])
    write(os.path.join(packet, 'submission', 'method.md'), '# Method\n')
    write(os.path.join(packet, 'new_file.txt'), 'x\n')
    code, out, g, r, _ = extract(tmp, packet, 'e_new', reads + [STEP, answer(GOOD)])
    check('grade: a new file in the packet rejects the grade', code != 0 and g is None, out[-300:])


if __name__ == '__main__':
    with tempfile.TemporaryDirectory() as tmp:
        packet = packet_tests(tmp)
        if packet is None:                                   # the packet tool failed: build a plain packet for the grade tests
            plain = os.path.join(tmp, 'plain')
            write(os.path.join(plain, 'results.json'), '{}\n')
            write(os.path.join(plain, 'method.md'), '# Method\n')
            packet = os.path.join(tmp, 'P0')
            code, out = run('make_packet.py', plain, os.path.join(tmp, 'score'), packet)
            if code != 0:
                sys.exit('no packet could be made: ' + out)
        grade_tests(tmp, packet)
    bad = results.count(False)
    print('%d of %d tests failed' % (bad, len(results)) if bad else 'ALL %d TOOL TESTS PASSED' % len(results))
    sys.exit(1 if bad else 0)
