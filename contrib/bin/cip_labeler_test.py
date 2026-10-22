"""Regression tests for original-CIP isotope and element labelling."""
import contextlib
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock

from rdkit import Chem

# run_python.sh uses -P; load the adjacent script without relying on sys.path.
_SCRIPT = Path(__file__).with_name("cip_labeler.py")
_SPEC = importlib.util.spec_from_file_location("cip_labeler", _SCRIPT)
labeler = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = labeler
_SPEC.loader.exec_module(labeler)


class CipLabelerTest(unittest.TestCase):
    def modified_molecule(self, smiles, options):
        original_writer = Chem.MolToSmiles
        captured = []

        def write(mol):
            captured.append(Chem.Mol(mol))
            return original_writer(mol)

        with mock.patch.object(Chem, "MolToSmiles", side_effect=write):
            output = labeler.cip_label(options, smiles)
        self.assertIsNotNone(Chem.MolFromSmiles(output, sanitize=False), output)
        return captured[0]

    def test_default_isotopes_use_original_assignments(self):
        smiles = "C[C@H](O)[C@H](F)Cl"
        centres = labeler.find_cip_centres(Chem.MolFromSmiles(smiles))
        self.assertEqual({c[2] for c in centres}, {'R', 'S'})
        mol = self.modified_molecule(smiles, labeler.Options())
        for index, element, label in centres:
            self.assertEqual(mol.GetAtomWithIdx(index).GetAtomicNum(), element)
            self.assertEqual(mol.GetAtomWithIdx(index).GetIsotope(),
                             8 if label == 'R' else 9)

    def test_element_replacements_use_snapshot(self):
        smiles = "C[C@H](O)[C@H](F)Cl"
        original = Chem.MolFromSmiles(smiles)
        centres = labeler.find_cip_centres(original)
        options = labeler.Options(r_elements={6: 24}, s_elements={6: 55})
        with mock.patch.object(labeler.rdCIPLabeler, 'AssignCIPLabels',
                               wraps=labeler.rdCIPLabeler.AssignCIPLabels) as assign:
            mol = self.modified_molecule(smiles, options)
            self.assertEqual(assign.call_count, 1)
        for index, _, label in centres:
            atom = mol.GetAtomWithIdx(index)
            self.assertEqual(atom.GetAtomicNum(), 24 if label == 'R' else 55)
            self.assertEqual(atom.GetIsotope(), 0)
            self.assertEqual(atom.GetTotalNumHs(), original.GetAtomWithIdx(index).GetTotalNumHs())
        for atom in original.GetAtoms():
            if not any(atom.GetIdx() == c[0] for c in centres):
                self.assertEqual(mol.GetAtomWithIdx(atom.GetIdx()).GetAtomicNum(), atom.GetAtomicNum())

    def test_unmapped_centres_fall_back_to_isotopes(self):
        smiles = "C[C@H](O)[C@H](F)Cl"
        mol = self.modified_molecule(smiles, labeler.Options(r_elements={6: 24}))
        for index, _, label in labeler.find_cip_centres(Chem.MolFromSmiles(smiles)):
            atom = mol.GetAtomWithIdx(index)
            self.assertEqual(atom.GetAtomicNum(), 24 if label == 'R' else 6)
            self.assertEqual(atom.GetIsotope(), 0 if label == 'R' else 9)

    def test_no_centres_returns_input_unchanged(self):
        options = labeler.Options(r_elements={6: 24})
        for smiles in ("OCC", "CC(F)Cl", "C/C=C/C"):
            self.assertEqual(labeler.cip_label(options, smiles), smiles)
        self.assertEqual(options.molecules_read, 3)
        self.assertEqual(options.molecules_with_chiral_centres, 0)

    def test_existing_isotope_survives_element_replacement(self):
        smiles = "C[13C@H](F)Cl"
        options = labeler.Options(r_elements={6: 24}, s_elements={6: 55})
        mol = self.modified_molecule(smiles, options)
        self.assertEqual(mol.GetAtomWithIdx(1).GetIsotope(), 13)

    def test_mapping_validation(self):
        self.assertEqual(labeler.parse_element_mappings(['C:Cr', 'N:Nd', 'C:Cr']),
                         {6: 24, 7: 60})
        for values in (['C'], ['C:'], ['C:Cr:Cs'], ['Xx:Cr'], ['C:*'], ['C:Cr', 'C:Cs']):
            with self.subTest(values=values), self.assertRaises(ValueError):
                labeler.parse_element_mappings(values)

    def test_cli_smiles_filter_and_descriptors(self):
        smiles = "C[C@H](O)[C@H](F)Cl"
        args = [sys.executable, str(_SCRIPT), '--R_element=C:Cr', '--S_element=C:Cs']
        output = subprocess.run([*args, '-'], input=smiles + ' sample name\n',
                                text=True, capture_output=True, check=True).stdout
        self.assertTrue(output.rstrip().endswith('sample name'))
        self.assertIn('Cr', output)
        self.assertIn('Cs', output)
        tdt = f'$SMI<{smiles}>\nPCN<sample>\n|\n'
        filtered = subprocess.run([*args, '--f', '-'], input=tdt, text=True,
                                  capture_output=True, check=True).stdout
        self.assertIn('Cr', filtered)
        self.assertTrue(filtered.endswith('PCN<sample>\n|\n'))
        options = labeler.Options(r_elements={6: 24}, s_elements={6: 55})
        with contextlib.redirect_stdout(io.StringIO()) as output:
            labeler.cip_label_descriptors(options, smiles + ' sample')
        self.assertEqual(output.getvalue(), 'sample 1 1\n')
        bad = subprocess.run([sys.executable, str(_SCRIPT), '--R_element=C:Xx', '-'],
                             input=smiles, text=True, capture_output=True)
        self.assertNotEqual(bad.returncode, 0)


if __name__ == '__main__':
    unittest.main()
