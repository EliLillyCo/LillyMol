"""Run with the current nanobind module on PYTHONPATH."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import medchem_reaction_overlap as overlap
import lillymol

ROOT = Path(__file__).resolve().parents[3]


class OverlapTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def write(self, name, text):
        path = self.directory / name
        path.write_text(text)
        return path

    def reaction(self, name, smarts='[C]-[OD1]'):
        return self.write(name, f'''name: "{name}"
scaffold {{ smarts: "{smarts}" change_element {{ atom: 1 element: "N" }} }}
''')

    def test_identity_containment_partial_and_no_overlap(self):
        records = overlap.read_molecules([self.write('in.smi', 'CO same id\nCCO same id\n')])
        a, b, c = (0, 'CN'), (1, 'CCN'), (1, 'CCC')
        for candidate, existing, expected in [
            ({a, b}, {a, b}, 'identical_on_corpus'),
            ({a}, {a, b}, 'candidate_contained_on_corpus'),
            ({a, b}, {a}, 'existing_contained_on_corpus'),
            ({a, b}, {a, c}, 'partial_overlap'),
            ({a}, {b}, 'no_overlap')]:
            with self.subTest(expected=expected):
                row = overlap.compare(candidate, existing, records, 1)
                self.assertEqual(row['relation'], expected)
                self.assertTrue(all(len(v) <= 1 for v in row['examples'].values()))
        self.assertEqual(overlap.compare({a}, {b}, records, 1)['shared'], 0)

    def test_duplicate_names_multiple_sites_and_self_skip(self):
        candidate = self.reaction('candidate.rxn')
        self.reaction('copy.rxn')
        library = self.write('REACTIONS', 'PROTO:candidate.rxn\nPROTO:copy.rxn\nPROTO:copy.rxn\n')
        inputs = self.write('in.smi', 'CO duplicate name\nCO duplicate name\nOCCO multisite\n')
        report = overlap.check(candidate, library, [inputs])
        self.assertEqual(report['candidate_products'], 3)
        self.assertEqual(report['reactions_checked'], 1)
        self.assertEqual(report['results'][0]['relation'], 'identical_on_corpus')
        self.assertEqual(report['errors'], [])
        mol = lillymol.MolFromSmiles('OCCO multisite')
        before = mol.unique_smiles()
        products = overlap.load_reaction(candidate, True).process(mol)
        self.assertEqual(len(products), 2)  # Both sites processed before deduplication.
        self.assertEqual(mol.unique_smiles(), before)

    def test_relative_query_file_and_legacy_loader(self):
        self.write('query.textproto', 'query { smarts: "[CD1]-[OD1]" }\n')
        path = self.write('relative.rxn', '''name: "relative"
scaffold { query_file: "query.textproto" change_element { atom: 1 element: "N" } }
''')
        wizard = overlap.load_reaction(path, True)
        listed = lillymol.MedchemWizard()
        listed.read_reactions(str(self.write('list', 'PROTO:relative.rxn\n')))
        mol = lillymol.MolFromSmiles('CO methanol')
        self.assertEqual(len(listed.process(mol)), 1)
        self.assertEqual(wizard.process(mol)[0].unique_smiles(),
                         lillymol.MolFromSmiles('CN').unique_smiles())
        reaction = lillymol.Reaction()
        self.assertTrue(reaction.construct_from_smirks('[C:1][O:2]>>[C:1][N:2]'))
        msi = self.directory / 'legacy.rxn'
        self.assertTrue(reaction.write_msi(str(msi)))
        legacy = overlap.load_reaction(msi, False)
        self.assertEqual(legacy.process(mol)[0].unique_smiles(),
                         wizard.process(mol)[0].unique_smiles())
        # Appending a bad file leaves loaded reactions and hit counters usable.
        with self.assertRaises(RuntimeError):
            legacy.read_reaction(str(self.write('bad.rxn', 'not a reaction')))
        self.assertEqual(legacy.number_reactions(), 1)
        self.assertEqual(len(legacy.process(mol)), 1)
        legacy.read_reaction(str(path))
        self.assertEqual(legacy.number_reactions(), 2)
        self.assertEqual(len(legacy.process(mol)), 2)

    def test_no_candidate_products_invalid_input_and_incomplete_library(self):
        candidate = self.reaction('candidate.rxn')
        library = self.write('REACTIONS', 'PROTO:missing.rxn\n')
        inputs = self.write('in.smi', 'CO id\n')
        report = overlap.check(candidate, library, [inputs])
        self.assertEqual(len(report['errors']), 1)
        self.assertEqual(report['reactions_checked'], 0)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(overlap.main([str(candidate), str(inputs), '--library', str(library)]), 2)
        inputs.write_text('CC id\n')
        with self.assertRaisesRegex(ValueError, 'no products'):
            overlap.check(candidate, library, [inputs])
        inputs.write_text('invalidSMILES id\n')
        with self.assertRaisesRegex(ValueError, ':1: invalid SMILES'):
            overlap.read_molecules([inputs])

    def test_real_amide_partial_overlap(self):
        directory = ROOT / 'data/MedchemWizard'
        candidate = directory / 'amide_to_1_2_4_oxadiazole.rxn'
        old = directory / 'Ester-amide_isosteres.1-3-4_oxadiazole.textproto'
        library = self.write('REACTIONS', f'PROTO:{old}\n')
        inputs = self.write('in.smi', 'O=C(NC)c1ccccc1 amide\nCC(=O)OC ester\n')
        report = overlap.check(candidate, library, [inputs])
        row = report['results'][0]
        self.assertEqual(row['shared'], 1)
        self.assertEqual(row['existing_only'], 1)
        self.assertEqual(row['relation'], 'candidate_contained_on_corpus')

    def test_json_examples_and_output_protection(self):
        candidate = self.reaction('candidate.rxn')
        self.reaction('copy.rxn')
        library = self.write('REACTIONS', 'PROTO:copy.rxn\n')
        inputs = self.write('in.smi', 'CO one\nCO two\n')
        output = self.directory / 'report.json'
        with contextlib.redirect_stdout(io.StringIO()):
            status = overlap.main([str(candidate), str(inputs), '--library',
                                   str(library), '--json', str(output),
                                   '--max-examples', '1'])
        self.assertEqual(status, 0)
        report = json.loads(output.read_text())
        self.assertEqual(report['results'][0]['shared'], 2)
        self.assertEqual(len(report['results'][0]['examples']['shared']), 1)
        before = candidate.read_text()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            overlap.main([str(candidate), str(inputs), '--library', str(library),
                          '--json', str(candidate)])
        self.assertEqual(candidate.read_text(), before)

    def test_binding_products_agree_with_trxn(self):
        executable = ROOT / 'src/bazel-bin/Molecule_Tools/trxn'
        if not executable.exists():
            self.skipTest('trxn not built')
        reaction = ROOT / 'data/MedchemWizard/amide_to_1_2_4_oxadiazole.rxn'
        inputs = ROOT / 'test/trxn/amide_to_1_2_4_oxadiazole/in/input.smi'
        subprocess.run([str(executable), '-P', str(reaction), '-m', 'each', '-z', 'i',
                        '-S', str(self.directory / 'products'), str(inputs)],
                       check=True, capture_output=True)
        cli = {lillymol.MolFromSmiles(line).unique_smiles()
               for line in (self.directory / 'products.smi').read_text().splitlines()}
        bound = {smiles for _, smiles in overlap.generate_products(
            overlap.load_reaction(reaction, True), overlap.read_molecules([inputs]))}
        self.assertEqual(cli, bound)


if __name__ == '__main__':
    unittest.main()
