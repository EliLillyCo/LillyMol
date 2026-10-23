import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from lillymol import *
from lillymol_tools import *


def _skip_guards(test_case):
    """Return True if the test should be skipped due to missing dependencies."""
    try:
        import numpy  # pylint: disable=import-outside-toplevel,unused-import
    except ImportError:
        test_case.skipTest("NumPy is not installed")
    if "LILLYMOL_HOME" not in os.environ:
        test_case.skipTest("LILLYMOL_HOME is not defined")


def _make_molecule(smiles, name):
    mol = Molecule()
    if not mol.build_from_smiles(f"{smiles} {name}"):
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


# ---------------------------------------------------------------------------
# Existing single-molecule tests (unchanged)
# ---------------------------------------------------------------------------

class TestDescriptorConvenienceFunctions(unittest.TestCase):

    def test_molecular_weight(self):
        mol = _make_molecule("CCO", "ethanol")
        self.assertAlmostEqual(molecular_weight(mol), 46.069, delta=0.001)

    def test_tpsa(self):
        mol = _make_molecule("CCO", "ethanol")
        value = tpsa(mol)
        self.assertIsNotNone(value)
        self.assertGreater(value, 0.0)

    def test_alogp_and_xlogp(self):
        mol = _make_molecule("CCO", "ethanol")
        self.assertIsNotNone(alogp(mol))
        mol = _make_molecule("CCO", "ethanol")
        self.assertIsNotNone(xlogp(mol))


class TestMolecularDescriptors(unittest.TestCase):

    def setUp(self):
        _skip_guards(self)
        self._descriptors = MolecularDescriptors()

    def test_names(self):
        names = self._descriptors.names()
        self.assertGreater(len(names), 0)
        self.assertEqual(names, self._descriptors.feature_names())
        self.assertEqual(names[0], "natoms")
        self.assertEqual(names[1], "nrings")

    def test_compute_array(self):
        mol = _make_molecule("CCO", "ethanol")
        values = self._descriptors.compute_array(mol)
        self.assertEqual(values.dtype.name, "float32")
        self.assertEqual(values.ndim, 1)
        self.assertEqual(values.size, len(self._descriptors.names()))

    def test_compute_dict(self):
        mol = _make_molecule("CCO", "ethanol")
        values = self._descriptors.compute(mol)
        self.assertIsInstance(values, dict)
        self.assertEqual(values["natoms"], 3.0)
        self.assertEqual(values["nrings"], 0.0)
        self.assertIn("amw", values)

    def test_compute_list(self):
        mols = [
            _make_molecule("CCO", "ethanol"),
            _make_molecule("c1ccccc1", "benzene"),
        ]
        values = self._descriptors.compute_list(mols)
        self.assertEqual(values.shape, (2, len(self._descriptors.names())))


class TestIWDescr(unittest.TestCase):

    def test_all_descriptors(self):
        try:
            import numpy  # pylint: disable=import-outside-toplevel,unused-import
        except ImportError:
            self.skipTest("NumPy is not installed")

        if "LILLYMOL_HOME" not in os.environ:
            self.skipTest("LILLYMOL_HOME is not defined")

        mol = Molecule()
        self.assertTrue(mol.build_from_smiles("CCO"))

        iwdescr = IWDescr()
        names = iwdescr.feature_names()
        values = iwdescr.process(mol)

        self.assertGreater(len(names), 0)
        self.assertEqual(len(names), len(set(names)))
        self.assertTrue(all(names))

        self.assertEqual(values.ndim, 1)
        self.assertEqual(values.dtype.name, "float32")
        self.assertEqual(values.size, len(names))

        self.assertEqual(names[0], "natoms")
        self.assertEqual(values[0], 3)   # natoms

        self.assertEqual(names[1], "nrings")
        self.assertEqual(values[1], 0)   # nrings


# ---------------------------------------------------------------------------
# process_list — NumPy array output
# ---------------------------------------------------------------------------

class TestProcessListNumpy(unittest.TestCase):

    def setUp(self):
        _skip_guards(self)
        self._iwdescr = IWDescr()
        self._molecules = [
            _make_molecule("CCO",        "ethanol"),
            _make_molecule("c1ccccc1",   "benzene"),
            _make_molecule("CC(=O)O",    "acetic_acid"),
            _make_molecule("c1ccccc1O",  "phenol"),
        ]

    def test_return_type_is_numpy_array(self):
        import numpy as np
        result = self._iwdescr.process_list(self._molecules)
        self.assertIsInstance(result, np.ndarray)

    def test_shape(self):
        result = self._iwdescr.process_list(self._molecules)
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.ndim, 2)
        self.assertEqual(result.shape, (len(self._molecules), ndescr))

    def test_dtype_is_float32(self):
        result = self._iwdescr.process_list(self._molecules)
        self.assertEqual(result.dtype.name, "float32")

    def test_row_matches_single_process(self):
        """Each row of process_list must equal the result of process() for that molecule."""
        import numpy as np
        # Re-build molecules so they have not been modified by setUp's process_list call.
        molecules = [
            _make_molecule("CCO",        "ethanol"),
            _make_molecule("c1ccccc1",   "benzene"),
            _make_molecule("CC(=O)O",    "acetic_acid"),
            _make_molecule("c1ccccc1O",  "phenol"),
        ]
        batch = self._iwdescr.process_list(molecules)

        for i, smi_name in enumerate([
            ("CCO",        "ethanol"),
            ("c1ccccc1",   "benzene"),
            ("CC(=O)O",    "acetic_acid"),
            ("c1ccccc1O",  "phenol"),
        ]):
            mol = _make_molecule(*smi_name)
            single = self._iwdescr.process(mol)
            np.testing.assert_array_equal(
                batch[i], single,
                err_msg=f"Row {i} ({smi_name[1]}) differs from process()"
            )

    def test_column_order_matches_feature_names(self):
        """Descriptor values are in the same column order as feature_names()."""
        names = self._iwdescr.feature_names()
        result = self._iwdescr.process_list(self._molecules)
        name_to_col = {n: c for c, n in enumerate(names)}

        # benzene (row 1) has 1 ring
        self.assertEqual(result[1, name_to_col["nrings"]], 1.0)
        # ethanol (row 0) has 0 rings
        self.assertEqual(result[0, name_to_col["nrings"]], 0.0)

    def test_empty_list_returns_zero_rows(self):
        """An empty molecule list should return a (0, ndescr) array, not an error."""
        result = self._iwdescr.process_list([])
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.shape, (0, ndescr))

    def test_single_molecule_list(self):
        """A one-element list should return a (1, ndescr) array."""
        mol = _make_molecule("CCO", "ethanol")
        result = self._iwdescr.process_list([mol])
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.shape, (1, ndescr))

    def test_default_does_not_return_dataframe(self):
        """Without as_dataframe=True the result must not be a DataFrame."""
        try:
            import pandas as pd  # pylint: disable=import-outside-toplevel
        except ImportError:
            return  # pandas not installed; nothing to check
        result = self._iwdescr.process_list(self._molecules)
        self.assertNotIsInstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# process_list — pandas DataFrame output
# ---------------------------------------------------------------------------

class TestProcessListDataFrame(unittest.TestCase):

    def setUp(self):
        _skip_guards(self)
        try:
            import pandas  # pylint: disable=import-outside-toplevel,unused-import
        except ImportError:
            self.skipTest("pandas is not installed")
        self._iwdescr = IWDescr()
        self._smiles_names = [
            ("CCO",        "ethanol"),
            ("c1ccccc1",   "benzene"),
            ("CC(=O)O",    "acetic_acid"),
            ("c1ccccc1O",  "phenol"),
        ]
        self._molecules = [_make_molecule(s, n) for s, n in self._smiles_names]

    def test_return_type_is_dataframe(self):
        import pandas as pd
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        self.assertIsInstance(result, pd.DataFrame)

    def test_shape(self):
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        ndescr = len(self._iwdescr.feature_names())
        self.assertEqual(result.shape, (len(self._molecules), ndescr))

    def test_column_names_match_feature_names(self):
        names = self._iwdescr.feature_names()
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        self.assertEqual(list(result.columns), names)

    def test_row_index_contains_molecule_names(self):
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        expected_index = [name for _, name in self._smiles_names]
        self.assertEqual(list(result.index), expected_index)

    def test_values_match_numpy_output(self):
        """DataFrame values must be numerically identical to the numpy array output."""
        import numpy as np
        # Rebuild so neither call's molecule-mutation affects the other.
        mols_np = [_make_molecule(s, n) for s, n in self._smiles_names]
        mols_df = [_make_molecule(s, n) for s, n in self._smiles_names]
        arr = self._iwdescr.process_list(mols_np)
        df  = self._iwdescr.process_list(mols_df, as_dataframe=True)
        np.testing.assert_array_equal(arr, df.values)

    def test_lookup_by_column_name(self):
        """DataFrame columns are accessible by descriptor name."""
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        # benzene has 1 ring; ethanol has 0
        self.assertEqual(result.loc["benzene",  "nrings"], 1.0)
        self.assertEqual(result.loc["ethanol",  "nrings"], 0.0)

    def test_lookup_by_row_index(self):
        """DataFrame rows are accessible by molecule name."""
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        row = result.loc["ethanol"]
        self.assertEqual(len(row), len(self._iwdescr.feature_names()))

    def test_dtype_is_float32(self):
        result = self._iwdescr.process_list(self._molecules, as_dataframe=True)
        self.assertEqual(result["natoms"].dtype.name, "float32")

    def test_empty_list_returns_empty_dataframe(self):
        """An empty molecule list should return a DataFrame with 0 rows and correct columns."""
        result = self._iwdescr.process_list([], as_dataframe=True)
        names = self._iwdescr.feature_names()
        self.assertEqual(result.shape[0], 0)
        self.assertEqual(list(result.columns), names)

    def test_single_molecule_list(self):
        """A one-element list should return a single-row DataFrame."""
        mol = _make_molecule("CCO", "ethanol")
        result = self._iwdescr.process_list([mol], as_dataframe=True)
        self.assertEqual(result.shape[0], 1)
        self.assertEqual(result.index[0], "ethanol")

    def test_false_flag_returns_numpy_not_dataframe(self):
        """Explicitly passing as_dataframe=False must return a numpy array."""
        import numpy as np
        import pandas as pd
        result = self._iwdescr.process_list(self._molecules, as_dataframe=False)
        self.assertIsInstance(result, np.ndarray)
        self.assertNotIsInstance(result, pd.DataFrame)


class TestJWCats(unittest.TestCase):

    def setUp(self):
        _skip_guards(self)

    def test_jwcats_default(self):
        import numpy as np
        jwcats = JWCats()
        names = jwcats.feature_names()
        self.assertEqual(len(names), 150)
        self.assertEqual(names[0], "jwc_B1PAA")
        self.assertEqual(names[-1], "jwc_B10PHH")

        mol = _make_molecule("CCO", "ethanol")
        values = jwcats.process(mol)
        self.assertIsInstance(values, np.ndarray)
        self.assertEqual(values.dtype.name, "float64")
        self.assertEqual(values.ndim, 1)
        self.assertEqual(values.size, len(names))

    def test_jwcats_without_hydrophobe_pairs(self):
        jwcats = JWCats()
        jwcats.set_include_hydrophobic_pairs(0)
        self.assertFalse(jwcats.initialised())
        self.assertTrue(jwcats.initialise())
        names = jwcats.feature_names()
        self.assertEqual(len(names), 140)
        self.assertNotIn("jwc_B1PHH", names)
        self.assertNotIn("jwc_B10PHH", names)

        mol = _make_molecule("CC", "ethane")
        values = jwcats.process(mol)
        self.assertEqual(values.size, len(names))

    def test_jwcats_process_list(self):
        jwcats = JWCats()
        mols = [
            _make_molecule("CCO", "ethanol"),
            _make_molecule("c1ccccc1", "benzene"),
        ]
        values = jwcats.process_list(mols)
        self.assertEqual(values.shape, (2, len(jwcats.feature_names())))

    def test_jwcats_can_be_built_without_default_assigners(self):
        jwcats = JWCats(False)
        self.assertTrue(jwcats.initialised())
        names = jwcats.feature_names()
        self.assertEqual(len(names), 150)

        mol = _make_molecule("CC", "ethane")
        values = jwcats.process(mol)
        name_to_col = {name: i for i, name in enumerate(names)}
        self.assertAlmostEqual(values[name_to_col["jwc_B1PHH"]], 0.5)


class TestSmallDescriptorHelpers(unittest.TestCase):

    def test_tpsa_compute(self):
        mol = _make_molecule("CCO", "ethanol")
        tpsa = TPSA()
        value = tpsa.compute(mol)
        self.assertIsNotNone(value)
        self.assertIsInstance(value, float)
        self.assertGreater(value, 0.0)

    def test_tpsa_empty_molecule_returns_none(self):
        tpsa = TPSA()
        self.assertIsNone(tpsa.compute(Molecule()))

    def test_tpsa_setters_are_available(self):
        tpsa = TPSA()
        tpsa.set_display_psa_unclassified_atom_messages(0)
        tpsa.set_return_zero_for_unclassified_atoms(0)
        tpsa.set_non_zero_contribution_for_SD2(1)
        tpsa.set_zero_for_all_sulphur_atoms(0)
        tpsa.set_zero_for_all_phosphorus_atoms(0)
        tpsa.set_convert_to_charge_separated(0)
        mol = _make_molecule("CCO", "ethanol")
        self.assertIsNotNone(tpsa.compute(mol))

    def test_hba_hbd(self):
        ethanol = _make_molecule("CCO", "ethanol")
        self.assertEqual(HbaHbd(ethanol), (1, 1))

        methylamine = _make_molecule("CN", "methylamine")
        self.assertEqual(HbaHbd(methylamine), (1, 2))

        acetic_acid = _make_molecule("CC(=O)O", "acetic_acid")
        self.assertEqual(HbaHbd(acetic_acid), (2, 1))


if __name__ == "__main__":
    unittest.main()
