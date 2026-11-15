"""Reference answers of the variant without the interface-crack definitions (task/VARIANT_NO_DEFINITIONS.md), written to
reference_no_definitions.json. They are the key's own numbers (reference.json), renamed, plus one new number:

- the amplitude A = lim sqrt(2 pi r) |T(r)| of the crack-tip traction T = sigma22 + i sigma12 is |K|;
- the phase of the crack-tip field at r = h is arg(K h^(i eps)) = psi, because T = K r^(i eps)/sqrt(2 pi r);
- when r increases tenfold, arg T changes by eps ln 10 (radians), here -8.514 degrees.
Teng, 14 Nov 2026: "How about we remove the interface-crack physics"; the wording of the variant approved on the same
day ("yes to all"). Run: python3 -I reference_no_definitions.py   (standard library)
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
REF = json.load(open(os.path.join(HERE, 'reference.json')))
out = {'variant': 'no_definitions', 'source': 'reference.json (phase angles from %s)' % REF['phase_angles_from'],
       'caseA': {'G_E2_h_over_P2': REF['caseA']['G_E2_h_over_P2'], 'amplitude_sqrt_h_over_P': REF['caseA']['Kabs_sqrt_h_over_P'],
                 'phase_at_h_deg': REF['caseA']['psi_deg']},
       'caseB': {'G_E2_h3_over_M2': REF['caseB']['G_E2_h3_over_M2'], 'amplitude_h32_over_M': REF['caseB']['Kabs_h32_over_M'],
                 'phase_at_h_deg': REF['caseB']['psi_deg']},
       'phase_at_h_both_deg': REF['psi_both_deg'],
       'phase_change_per_decade_deg': math.degrees(REF['eps'] * math.log(10.0))}
assert abs(out['phase_change_per_decade_deg'] + 8.5136) < 1e-3
json.dump(out, open(os.path.join(HERE, 'reference_no_definitions.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
