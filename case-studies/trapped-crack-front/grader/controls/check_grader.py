"""Compare the grader's answers on the control submissions with EXPECTED.json.

usage: python3 -I check_grader.py RUN_DIR      RUN_DIR holds MAPPING.json ({"G1": "c3_rotation_theory", ...}) and, for
each packet name, the folder <name>.grade/ written by extract_grade.py. Exit 1 if an answer differs or a grade is missing.
The first wrong step named for a wrong submission is not compared here: a person reads it.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
exp = json.load(open(os.path.join(HERE, 'EXPECTED.json')))
run = os.path.abspath(sys.argv[1])
bad = 0
for packet, control in sorted(json.load(open(os.path.join(run, 'MAPPING.json'))).items()):
    p = os.path.join(run, packet + '.grade', 'grade.json')
    if not os.path.isfile(p):
        print('FAIL %s (%s): no grade.json' % (packet, control)); bad += 1; continue
    g, e, diffs = json.load(open(p)), exp[control], []
    for key in ('theory_source', 'claimed_checks', 'error_location'):
        if g[key]['answer'] != e[key]:
            diffs.append('%s: %s (expected %s)' % (key, g[key]['answer'], e[key]))
    for c, want in e['coefficients_equal_key'].items():
        got = g['theory_source']['coefficients'][c]['equals_key']
        if got != want:
            diffs.append('coefficient %s equals_key: %s (expected %s)' % (c, got, want))
    n_uns = sum(c.get('status') == 'unsupported' for c in g['claimed_checks']['checks'])
    if n_uns < e.get('unsupported_at_least', 0):
        diffs.append('%d checks marked unsupported (expected at least %d)' % (n_uns, e['unsupported_at_least']))
    bad += bool(diffs)
    print('%s %s (%s): %s' % ('FAIL' if diffs else 'ok  ', packet, control, '; '.join(diffs) if diffs else
          '%s | %s | %s' % tuple(g[k]['answer'] for k in ('theory_source', 'claimed_checks', 'error_location'))))
    if g['error_location']['answer'] not in ('none',):
        print('       first wrong step: %s  [%s]' % (g['error_location'].get('first_wrong_step'), g['error_location'].get('where')))
print('FAILED %d' % bad if bad else 'GRADER AGREES WITH ALL CONTROLS')
sys.exit(1 if bad else 0)
