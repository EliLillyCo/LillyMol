"""Run with: ../../run_python.sh -m unittest discover -s test -p medchem_transform_activity_html_test.py"""
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import medchem_transform_activity_html as report


class ReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.starts = self.write('starts.smi', 'C p1\nCC p2\n')
        self.activity = self.write('activity.txt', 'id activity\np1 1\np2 2\nx 3\ny 0\n')
        self.found = self.write('found.smi', 'CCCC p2__MCTA__grow x\nCC p1__MCTA__grow x\nCO p1__MCTA__oxygen y\n')
        self.notfound = self.write('notfound.smi', 'CCC p2__MCTA__grow\nCN p1__MCTA__nitrogen\n')
        self.pred = self.write('notfound.pred', 'id prediction\np1__MCTA__nitrogen 4\np2__MCTA__grow 1\n')

    def write(self, name, value):
        path = self.root / name
        path.write_text(value)
        return path

    def build(self, **kwargs):
        return report.build_data(self.starts, self.activity, self.found, self.notfound, self.pred, **kwargs)

    def test_counts_changes_and_smallest_global_exemplar(self):
        data = self.build(max_examples=1)
        rows = {r['name']: r for r in data['rows']}
        grow = rows['grow']
        self.assertEqual(grow['observed']['n'], 2)
        self.assertEqual(grow['observed']['parents'], 2)
        self.assertEqual(grow['observed']['mean'], 1.5)
        self.assertEqual(grow['predicted']['median'], -1)
        self.assertEqual(grow['examples'][0]['after'], 'CC')
        self.assertEqual(grow['examples'][0]['before'], 'C')
        self.assertEqual(grow['examples'][0]['atoms'], 2)
        self.assertEqual(rows['oxygen']['observed']['median'], -1)
        self.assertEqual(rows['nitrogen']['predicted']['median'], 3)
        self.assertTrue(all(len(r['examples']) <= 1 for r in data['rows']))
        self.assertTrue(all('<svg' in svg for svg in data['drawings'].values()))

    def test_many_records_keep_all_statistics_but_cap_examples(self):
        self.pred.write_text('id prediction\n' + 'p1__MCTA__nitrogen 4\n' * 1000 + 'p2__MCTA__grow 1\n')
        data = self.build(max_examples=1)
        nitrogen = next(r for r in data['rows'] if r['name'] == 'nitrogen')
        self.assertEqual(nitrogen['predicted']['n'], 1000)
        self.assertEqual(nitrogen['predicted']['parents'], 1)
        self.assertEqual(sum(nitrogen['predicted']['hist']), 1000)
        self.assertEqual(nitrogen['predicted']['dots'], [])
        self.assertEqual(len(nitrogen['examples']), 1)

    def test_smallest_example_can_come_from_predictions(self):
        self.notfound.write_text('N p2__MCTA__grow\nCN p1__MCTA__nitrogen\n')
        grow = next(r for r in self.build(max_examples=1)['rows'] if r['name'] == 'grow')
        self.assertEqual(grow['examples'][0]['source'], 'predicted')
        self.assertEqual(grow['examples'][0]['atoms'], 1)

    def test_sorted_fast_path(self):
        self.found.write_text('CC p1__MCTA__grow x\nCCCC p2__MCTA__grow x\nCO p1__MCTA__oxygen y\n')
        a = self.build(max_examples=1)
        b = self.build(max_examples=1, products_sorted=True)
        self.assertEqual(a, b)

    def test_ambiguous_prediction_ids_and_equivalent_structures(self):
        self.notfound.write_text('CCC p2__MCTA__grow\nCO p1__MCTA__nitrogen\nCN p1__MCTA__nitrogen\n')
        nitrogen = next(r for r in self.build()['rows'] if r['name'] == 'nitrogen')
        self.assertEqual(nitrogen['predicted']['n'], 1)
        self.assertTrue(all(e['ambiguous'] and e['delta'] is None for e in nitrogen['examples']))
        self.notfound.write_text('CCC p2__MCTA__grow\nCO p1__MCTA__nitrogen\nOC p1__MCTA__nitrogen\n')
        nitrogen = next(r for r in self.build()['rows'] if r['name'] == 'nitrogen')
        self.assertEqual(len(nitrogen['examples']), 1)
        self.assertFalse(nitrogen['examples'][0]['ambiguous'])
        self.assertEqual(nitrogen['examples'][0]['delta'], 3)

    def test_legacy_inputs_and_csv_activity(self):
        self.found.write_text('CC p1 grow x\n')
        self.notfound.write_text('CN p1 nitrogen\n')
        self.pred.write_text('id transform prediction\np1 nitrogen 4\n')
        self.activity.write_text('id,activity\np1,1\np2,2\nx,3\ny,0\n')
        renamed = self.activity.with_suffix('.csv')
        self.activity.rename(renamed)
        self.activity = renamed
        self.assertEqual(len(self.build()['rows']), 2)

    def test_neutral_interval_and_dense_histogram(self):
        summary = report.summarize([-.2, 0, .1, 2], {'p1'}, 2, .2)
        self.assertEqual(summary['positive'], .25)
        self.assertEqual(summary['negative'], 0)
        dense = report.summarize(list(range(200)), {'p1', 'p2'}, 199, 0)
        self.assertEqual(sum(dense['hist']), 200)
        self.assertEqual(dense['dots'], [])

    def test_unscored_products_and_unmatched_predictions(self):
        self.notfound.write_text(self.notfound.read_text() + 'CC p1__MCTA__unscored\n')
        data = self.build()
        self.assertEqual(data['unscored_products'], 1)
        self.assertNotIn('unscored', [r['name'] for r in data['rows']])
        self.pred.write_text('id prediction\np1__MCTA__missing 4\n')
        with self.assertRaisesRegex(ValueError, 'no notfound structure'):
            self.build()

    def test_validation_and_lillymol_analysis_boundary(self):
        with mock.patch.object(report.lillymol, 'MolFromSmiles', wraps=report.lillymol.MolFromSmiles) as parse:
            report.structure_info.cache_clear()
            self.assertEqual(report.structure_info('C[C@H](O)N')[0], 4)
            self.assertTrue(parse.called)
        self.pred.write_text('id prediction\np1__MCTA__nitrogen nan\n')
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            self.build()

    def test_html_is_offline_and_escapes_input(self):
        data = self.build()
        data['rows'][0]['name'] = '</script><script>alert(1)</script>'
        page = report.render_html(data, '<unsafe>')
        self.assertIn('&lt;unsafe&gt;', page)
        self.assertNotIn('<script>alert(1)</script>', page)
        payload = re.search(r'<script id="data" type="application/json">(.*?)</script>', page, re.S).group(1)
        self.assertEqual(json.loads(payload)['rows'][0]['name'], data['rows'][0]['name'])
        self.assertNotRegex(page, r'<script[^>]+src=|<link[^>]+href=')

    def test_cli_and_input_overwrite_protection(self):
        output = self.root / 'report.html'
        args = ['--starting-smiles', str(self.starts), '--activity', str(self.activity),
                '--found', str(self.found), '-o', str(output)]
        self.assertIsNone(report.main(args))
        self.assertIn('Medchem transformation activity', output.read_text())
        if shutil.which('node'):
            harness = Path(__file__).with_name('medchem_transform_activity_html_ui_test.js')
            result = subprocess.run(['node', str(harness), str(output)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        original = self.found.read_bytes()
        args[-1] = str(self.found)
        with self.assertRaises(SystemExit):
            report.main(args)
        self.assertEqual(self.found.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
