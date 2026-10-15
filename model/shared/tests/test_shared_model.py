"""Focused tests for renamed kinetics, voltage timing, and condition resolution."""
import copy
import importlib.util
import json
from pathlib import Path
import re
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build import COMMON, MODEL_ROOT, load_condition, materialize, render_bngl, validate
from batch import suite_leaves


class SharedModelTests(unittest.TestCase):
    def test_all_conditions_and_suite_reuse(self):
        leaves = suite_leaves('mixed_distance') + suite_leaves('spont_distance')
        self.assertEqual(len(set(leaves)), 12)
        self.assertIn(suite_leaves('spont_1_vs_2')[1], suite_leaves('spont_distance'))
        for leaf in leaves:
            c = load_condition(leaf)
            self.assertEqual(c['duration_s'], .045)
            self.assertEqual(c['waveform_start_s'], .01)
            self.assertEqual(sum(x['glutamate'] for x in c['clusters']), 500)
            self.assertEqual([x['release_time_s'] for x in c['clusters']], [.01, .01])
            self.assertEqual(sum(c['vgcc_populations'].values()), 1)

    def test_nmda_mass_action_rates_and_population_identity(self):
        c = load_condition(suite_leaves('mixed_distance')[0])
        text = render_bngl(c)
        kinetics = json.loads((COMMON / 'nmda_kinetics.json').read_text())
        for population in ('sNMDAR', 'eNMDAR'):
            expected = []
            for rule in kinetics['rules']:
                lhs, rhs = rule.split('<->')
                m = re.fullmatch(r'\s*(.*?)\s+(\S+),\s*(\S+)\s*', rhs)
                product, fwd, rev = m.groups()
                expected.extend([f'{lhs.strip()} -> {product} {fwd}', f'{product} -> {lhs.strip()} {rev}'])
            actual = [line.split(': ', 1)[1].replace(population, 'NMDAR')
                      for line in text.splitlines() if line.strip().startswith(population + '_gating_')]
            self.assertEqual(actual, expected)
        for line in text.splitlines():
            if '->' in line:
                self.assertFalse('sNMDAR' in line and 'eNMDAR' in line)
        self.assertNotRegex(text, r'(?m)^\s*NMDAR(?:_[A-Za-z0-9]+)?\s*$')
        self.assertNotRegex(text, r'(?m)^\s*VSCC(?:_[0-3])?\s*$')

    def test_rate_tables_and_flux_rules(self):
        c = load_condition(suite_leaves('mixed_distance')[0])
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp) / 'snapshot'
            materialize(c, folder)
            fake_mcell = SimpleNamespace(Complex=lambda name, **kw: name,
                                         ReactionRule=lambda **kw: SimpleNamespace(**kw))
            spec = importlib.util.spec_from_file_location('tested_channels', folder / 'channels.py')
            channels = importlib.util.module_from_spec(spec)
            with patch.dict(sys.modules, mcell=fake_mcell,
                            condition=SimpleNamespace(CONFIG=c, PROTOCOL=c['channel_protocol'])):
                spec.loader.exec_module(channels)
                rules = []
                nmda, vgcc = channels.add_channel_rules(SimpleNamespace(add_reaction_rule=rules.append))
            self.assertEqual(len(rules), 15)  # 2 * (block, unblock, flux) + 9 VGCC rules
            self.assertEqual(set(nmda), {'sNMDAR', 'eNMDAR'})
            self.assertEqual(set(vgcc), {'eVGCC'})
            self.assertEqual(nmda['sNMDAR'].fwd_rate, channels.REST['k_NMDAR_Ca'])
            self.assertEqual(nmda['eNMDAR'].products, ['eNMDAR', 'Ca2@IN'])
            self.assertEqual(nmda['sNMDAR'].products, ['sNMDAR', 'Ca2@IN'])
            rates = nmda['eNMDAR'].variable_rate
            self.assertEqual(rates[0], [0., channels.REST['k_NMDAR_Ca']])
            for actual, original in zip(rates[1:], channels.EVOKED['nmdar_rate_scaled']):
                self.assertAlmostEqual(actual[0], original[0] + .01, places=12)
                self.assertEqual(actual[1], original[1])
            self.assertTrue(all(b[0] > a[0] for a,b in zip(rates,rates[1:])))
            for r in rules:
                self.assertFalse('eNMDAR' in r.reactants[0] and any('sNMDAR' in p for p in r.products))

    def test_invalid_layout_and_vgcc_duplication_rejected(self):
        config = load_condition(suite_leaves('mixed_distance')[0])
        c = copy.deepcopy(config)
        c['vgcc_populations']['eVGCC'] = 2
        with self.assertRaises(ValueError): validate(c)
        c = copy.deepcopy(config)
        c['clusters'][1]['walls'] = c['clusters'][0]['walls']
        with self.assertRaises(ValueError): validate(c)
        c = copy.deepcopy(config)
        c['clusters'][0]['population'] = 'NMDAR'
        with self.assertRaises(ValueError): validate(c)

    def test_snapshot_and_smoke_are_explicit(self):
        leaf = suite_leaves('mixed_distance')[0]
        original = load_condition(leaf)
        c = load_condition(leaf, smoke=True)
        self.assertEqual(original['duration_s'], .045)
        self.assertEqual(c['duration_s'], .0002)
        self.assertEqual(c['production_timeline']['duration_s'], .045)
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / 'snapshot'
            hashes = materialize(c, out)
            self.assertIn('geometry/Mushroom_dendrite.py', hashes)
            for path in out.rglob('*.py'):
                compile(path.read_text(), str(path), 'exec')
            with self.assertRaises(FileExistsError): materialize(c, out)


if __name__ == '__main__':
    unittest.main()
