"""Check kinetic equivalence, source conservation and absence of ER machinery."""
import json
from pathlib import Path
import re
import unittest
from build import ROOT, render, validate
from calcium_labels import CARRIERS, labeled, relabel
from geometry_checks import check_placement


def rules(text):
    block = text.split('begin reaction rules')[1].split('end reaction rules')[0]
    return dict(line.strip().split(': ', 1) for line in block.splitlines() if line.strip())


def inventory(side, carriers):
    return sum(len(re.findall(r'\b' + re.escape(c) + r'\b', side)) for c in carriers)


class ChemistryTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((ROOT / 'condition.json').read_text())
        self.text = render(self.config)
        self.rules = rules(self.text)

    def test_config_and_no_er(self):
        validate(self.config)
        active = self.text.split('begin molecule types')[1]
        self.assertIsNone(re.search(r'\b(?:SERCA\w*|RyR\w*|IP3R|Ca2_ER)\b', active))
        self.assertFalse((ROOT / 'geometry/SA.py').exists())
        self.assertNotIn('SA', (ROOT / 'runtime/instantiation.py').read_text())

    def test_every_label_rule_has_identical_unlabeled_rule(self):
        for name, rule in self.rules.items():
            if name.endswith(('_A', '_B')):
                tag = name[-1]
                self.assertEqual(rule, relabel(self.rules[name[:-2]], tag))
        for name, rule in self.rules.items():
            if name.startswith('NMDAR_') or name.endswith(('_A', '_B')):
                continue
            if any(re.search(r'\b' + c + r'\b', rule) for c in CARRIERS):
                if name in ('rxn_PMCA_leak', 'rxn_NCX_leak'):
                    self.assertNotIn(name + '_A', self.rules)
                else:
                    for tag in 'AB':
                        self.assertEqual(self.rules[name + '_' + tag], relabel(rule, tag))

    def test_label_conservation_except_explicit_sinks(self):
        for name, rule in self.rules.items():
            left, right = rule.split('->')
            for tag in 'AB':
                carriers = [labeled(c, tag) for c in CARRIERS]
                delta = inventory(right, carriers) - inventory(left, carriers)
                expected = -1 if name in [s + '_' + tag for s in ('Ca2Null','rxn_PMCA','rxn_NCX')] else 0
                self.assertEqual(delta, expected, name)

    def test_capacity_conservation(self):
        groups = [['Bm','CaBm'], ['Bf','CaBf'], ['PMCA','PMCA1'], ['NCX','NCX1'], ['PLC','PLC_Ca','PLC_Ca_PIP2']]
        for group in groups:
            expanded = group + [labeled(c,t) for c in group if c in CARRIERS for t in 'AB']
            for name, rule in self.rules.items():
                left, right = rule.split('->')
                self.assertEqual(inventory(left,expanded),inventory(right,expanded),name)

    def test_nmda_gating_identical_after_population_rename(self):
        for name, rule in self.rules.items():
            if name.startswith('NMDAR_A_gating'):
                self.assertEqual(rule.replace('NMDAR_A','NMDAR_B'), self.rules[name.replace('NMDAR_A','NMDAR_B')])

    def test_release_sites_are_in_cleft_for_both_distances(self):
        check_placement(self.config, ROOT/'geometry')
        config200 = json.loads((ROOT/'condition_200nm.json').read_text())
        placement = check_placement(config200, ROOT/'geometry')
        self.assertAlmostEqual(placement['center_distance_xy_nm'],200.083,places=3)
        config200['clusters'][1].pop('glutamate_z_um')
        with self.assertRaisesRegex(ValueError,'outside the local cleft'):
            check_placement(config200, ROOT/'geometry')


if __name__ == '__main__':
    unittest.main()
