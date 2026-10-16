import atexit
import copy
import math
import os
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from google.protobuf import text_format

import lillymol
from Molecule_Lib import substructure_pb2
from Utilities.GFP_Tools import nearneighbours_pb2


def _trace_process_exit():
    if os.environ.get("LILLYMOL_NB_TEST_TRACE"):
        print("ATEXIT unittest complete", file=sys.stderr, flush=True)


atexit.register(_trace_process_exit)


def _query_dir(kind):
    for envvar in ("C3TK_DATA_PERSISTENT", "LILLYMOL_HOME"):
        home = os.environ.get(envvar)
        if not home:
            continue
        if envvar == "C3TK_DATA_PERSISTENT":
            candidate = os.path.join(home, "queries", kind)
        else:
            candidate = os.path.join(home, "data", "queries", kind)
        if os.path.isdir(candidate):
            return candidate

    # Direct execution from the source checkout has nanobind/ two levels below
    # the repository root. Bazel runfiles only see declared data, so this path
    # normally exists only outside sandboxed test execution.
    candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                             "data", "queries", kind))
    if os.path.isdir(candidate):
        return candidate

    return None


def _hbonds_query_dir():
    return _query_dir("hbonds")


def _charges_query_dir():
    return _query_dir("charges")


def _qed_query_dir():
    return _query_dir("QED")


def _ring_replacement_file(fname="6a.smi"):
    home = os.environ.get("LILLYMOL_HOME")
    if home:
        candidate = os.path.join(home, "contrib", "test", "ring_replacement", fname)
        if os.path.exists(candidate):
            return candidate

    candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                             "contrib", "test", "ring_replacement", fname))
    if os.path.exists(candidate):
        return candidate

    return None


def _medchemwizard_reactions_file():
    candidates = []

    home = os.environ.get("LILLYMOL_HOME")
    if home:
        candidates.append(os.path.join(home, "data", "MedchemWizard", "REACTIONS"))

    candidates.extend([
        os.path.join(os.getcwd(), "data", "MedchemWizard", "REACTIONS"),
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                     "data", "MedchemWizard", "REACTIONS")),
    ])

    for envvar in ("TEST_SRCDIR", "RUNFILES_DIR"):
        root = os.environ.get(envvar)
        if root:
            candidates.extend([
                os.path.join(root, "_main", "data", "MedchemWizard", "REACTIONS"),
                os.path.join(root, "data", "MedchemWizard", "REACTIONS"),
            ])

    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate

    return None


def _methane():
    mol = lillymol.Molecule()
    assert mol.build_from_smiles("C methane")
    return mol


_CRC32C_TABLE = None


def _crc32c_table():
    global _CRC32C_TABLE
    if _CRC32C_TABLE is not None:
        return _CRC32C_TABLE
    table = []
    for i in range(256):
        crc = i
        for _ in range(8):
            if crc & 1:
                crc = (crc >> 1) ^ 0x82F63B78
            else:
                crc >>= 1
        table.append(crc & 0xffffffff)
    _CRC32C_TABLE = table
    return table


def _crc32c(data):
    crc = 0xffffffff
    table = _crc32c_table()
    for byte in data:
        crc = table[(crc ^ byte) & 0xff] ^ (crc >> 8)
    return (~crc) & 0xffffffff


def _masked_crc32c(data):
    crc = _crc32c(data)
    return (((crc >> 15) | ((crc << 17) & 0xffffffff)) + 0xa282ead8) & 0xffffffff


def _write_tfdatarecord(protos):
    fd, fname = tempfile.mkstemp(suffix=".nn.tfrecord")
    with os.fdopen(fd, "wb") as output:
        for proto in protos:
            data = proto.SerializeToString()
            length = struct.pack("<Q", len(data))
            output.write(length)
            output.write(struct.pack("<I", _masked_crc32c(length)))
            output.write(data)
            output.write(struct.pack("<I", _masked_crc32c(data)))
    return fname


def _nearneighbours(name, nbrs):
    proto = nearneighbours_pb2.NearNeighbours(name=name)
    for nbr_id, dist in nbrs:
        nbr = proto.nbr.add()
        nbr.id = nbr_id
        nbr.dist = dist
    return proto


def _nearneighbours_indices(name, nbrs):
    proto = nearneighbours_pb2.NearNeighboursIndices(name=name)
    for nbr_id, dist in nbrs:
        nbr = proto.nbr.add()
        nbr.id = nbr_id
        nbr.dist = dist
    return proto


def _write_text_file(contents, suffix=".txt"):
    fd, fname = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w") as output:
        output.write(contents)
    return fname


def _testdata_file(fname):
    candidates = []
    for envvar in ("TEST_SRCDIR", "RUNFILES_DIR"):
        root = os.environ.get(envvar)
        if root:
            candidates.extend([
                os.path.join(root, "_main", "nanobind", "testdata", fname),
                os.path.join(root, "nanobind", "testdata", fname),
            ])
    candidates.extend([
        os.path.join(os.path.dirname(__file__), "testdata", fname),
        os.path.join(os.path.dirname(__file__), "..", "nanobind", "testdata", fname),
        os.path.join(os.path.dirname(__file__), "..", "..", "nanobind", "testdata", fname),
    ])
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(fname)


class LillyMolNanobindTestCase(unittest.TestCase):

    def setUp(self):
        if os.environ.get("LILLYMOL_NB_TEST_TRACE"):
            print(f"START {self.id()}", file=sys.stderr, flush=True)

    def tearDown(self):
        if os.environ.get("LILLYMOL_NB_TEST_TRACE"):
            print(f"END {self.id()}", file=sys.stderr, flush=True)

CARBON_GFP = """
$SMI<C>
PCN<Methane>
FCTS<.E..........2;1;1;1;1>
|
$SMI<CC>
PCN<Ethane>
FCTS<.U..........2;1;1;1;1>
|
"""


SMILES = """C(=O)NCCCC CHEMBL45466
C(O)(=O)CCC(O)=O CHEMBL1200345
N1=NC(=CS1)C(O)=O CHEMBL247337
NOC1CCC(N)CC1 CHEMBL1213430
N1(=C(C)CCC1(C)C)=O CHEMBL325242
O=C(NC)[C@H]1NC(=O)CC1 CHEMBL1892080
C1=CC=C2C(=CC=N2)N1C CHEMBL593929
C1=CC=C2C(=C1)NC(N)S2 CHEMBL568765
C1C(N)CC1(O)P(O)(C)=O CHEMBL509338
N1(C(C#N)C1)C(=O)NCC CHEMBL150159
"""


class TestNanobindMolecule(LillyMolNanobindTestCase):

    def test_build_from_smiles(self):
        mol = lillymol.Molecule()
        self.assertTrue(mol.build_from_smiles("CCO ethanol"))
        self.assertEqual(mol.natoms(), 3)
        self.assertEqual(mol.nedges(), 2)
        self.assertEqual(mol.name(), "ethanol")
        self.assertEqual(mol.smiles(), "CCO")
        self.assertEqual(mol.unique_smiles(), "OCC")
        self.assertEqual(mol.molecular_formula(), "C2OH6")

    def test_mol_from_smiles(self):
        mol = lillymol.MolFromSmiles("c1ccccc1 benzene")
        self.assertIsNotNone(mol)
        self.assertEqual(mol.natoms(), 6)
        self.assertEqual(mol.nrings(), 1)
        self.assertEqual(mol.name(), "benzene")
        self.assertIsNone(lillymol.MolFromSmiles("["))

    def test_lillymol_from_smiles_and_batch_smiles(self):
        mol = lillymol.LillyMolFromSmiles("CCO ethanol")
        self.assertIsNotNone(mol)
        self.assertEqual(mol.name(), "ethanol")

        molecules = lillymol.MolFromSmiles(["C methane", "CC ethane", "["])
        self.assertEqual(len(molecules), 3)
        self.assertEqual(molecules[0].name(), "methane")
        self.assertEqual(molecules[1].natoms(), 2)
        self.assertEqual(molecules[2].natoms(), 0)

    def test_molecule_equality(self):
        mol1 = lillymol.MolFromSmiles("CCO ethanol")
        mol2 = lillymol.MolFromSmiles("CCO other")
        mol3 = lillymol.MolFromSmiles("CCN ethylamine")
        self.assertEqual(mol1, mol2)
        self.assertNotEqual(mol1, mol3)
        mol2.set_atomic_number(2, 7)
        self.assertEqual(mol2, mol3)

    def test_reader_open_and_next(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            # Readers must be closed before the temporary directory goes away.
            # On a local filesystem unlinking a file that is still open works, so
            # leaving it open passes; on NFS the unlink becomes a silly rename to
            # .nfsXXXX, the directory is then not empty and TemporaryDirectory
            # cleanup fails with ENOTEMPTY. The bazel sandbox puts /tmp on NFS
            # here, so a leaked reader fails under bazel and passes when run by
            # hand.
            reader = lillymol.Reader()
            self.assertTrue(reader.open(fname))
            mol = reader.next()
            self.assertIsNotNone(mol)
            self.assertEqual(mol.name(), "CHEMBL45466")
            self.assertEqual(reader.molecules_read(), 1)
            reader.close()

    def test_reader_bad_suffix_and_explicit_type(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.foo")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            reader = lillymol.Reader()
            self.assertFalse(reader.open(fname))
            self.assertTrue(reader.open(fname, lillymol.FileType.SMI))
            self.assertEqual(reader.next().name(), "CHEMBL45466")
            reader.close()

    def test_reader_iter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            reader = lillymol.Reader()
            self.assertTrue(reader.open(fname, lillymol.FileType.SMI))
            names = [mol.name() for mol in reader]
            self.assertEqual(len(names), 10)
            self.assertEqual(reader.molecules_read(), 10)
            reader.close()

    def test_slurp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            molecules = lillymol.slurp(fname)
            self.assertIsNotNone(molecules)
            self.assertEqual(len(molecules), 10)
            self.assertEqual(molecules[0].name(), "CHEMBL45466")
            self.assertEqual([mol.name() for mol in molecules[:2]], ["CHEMBL45466", "CHEMBL1200345"])

    def test_global_sdf_option_helpers_are_available(self):
        self.assertTrue(lillymol.set_sdf_identifier("idnumber"))
        lillymol.set_prepend_sdfid(True)
        lillymol.set_allsdfid(False)
        lillymol.set_sdf_tags_to_json(False)
        lillymol.set_firstsdftag(False)
        lillymol.set_ignore_bad_m(False)
        lillymol.set_mdlquiet(False)
        lillymol.set_allow_deuterium(False)
        lillymol.set_allow_tritium(False)

    def test_context_aliases(self):
        self.assertIs(lillymol.ReaderContext, lillymol.MolReaderContext)
        self.assertIs(lillymol.ContextWriter, lillymol.MolWriterContext)

    def test_mol_reader_context(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write(SMILES)

            with lillymol.MolReaderContext(fname, lillymol.FileType.SMI) as reader:
                names = [mol.name() for mol in reader]
            self.assertEqual(len(names), 10)
            self.assertEqual(names[0], "CHEMBL45466")

    def test_mol_reader_context_preprocessing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write("CC.O mixture\n")

            with lillymol.MolReaderContext(fname, largest_fragment=True) as reader:
                mol = reader.next()
                self.assertIsNotNone(mol)
                self.assertEqual(mol.natoms(), 2)
                self.assertIsNone(reader.next())

    def test_mol_reader_context_skip_invalid_valence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write("CC good\nC=[SH]-C bad\nCCC also_good\n")

            with lillymol.MolReaderContext(fname, skip_invalid_valence=True) as reader:
                self.assertEqual([mol.name() for mol in reader], ["good", "also_good"])
                self.assertEqual(reader.skipped_invalid_valence(), 1)
                self.assertEqual(reader.molecules_skipped(), 1)

    def test_mol_reader_context_skip_non_organic(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write("CC ethane\nCC[Na] sodium\nCCC propane\n")

            options = lillymol.ReaderOptions(skip_non_organic=True)
            with lillymol.MolReaderContext(fname, options=options) as reader:
                self.assertEqual([mol.name() for mol in reader], ["ethane", "propane"])
                self.assertEqual(reader.skipped_non_organic(), 1)
                self.assertEqual(reader.molecules_skipped(), 1)

    def test_mol_reader_context_preprocesses_before_skip_non_organic(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write("CC.[Na+] sodium_counterion\nCC[Na] covalent_sodium\n")

            with lillymol.MolReaderContext(fname, largest_fragment=True,
                                           skip_non_organic=True) as reader:
                molecules = list(reader)
                self.assertEqual([mol.name() for mol in molecules], ["sodium_counterion"])
                self.assertEqual(molecules[0].smiles(), "CC")
                self.assertEqual(reader.skipped_non_organic(), 1)
                self.assertEqual(reader.molecules_skipped(), 1)

    def test_mol_reader_context_skip_non_periodic_table_elements(self):
        lillymol.set_auto_create_new_elements(False)
        lillymol.set_atomic_symbols_can_have_arbitrary_length(False)
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.smi")
            with open(fname, "w") as writer:
                writer.write("CC ethane\nC[Zz] nonperiodic\nCCC propane\n")

            with lillymol.MolReaderContext(
                    fname, skip_non_periodic_table_elements=True) as reader:
                self.assertEqual([mol.name() for mol in reader], ["ethane", "propane"])
                self.assertEqual(reader.skipped_non_periodic_table_elements(), 1)
                self.assertEqual(reader.molecules_skipped(), 1)

        lillymol.set_auto_create_new_elements(False)
        lillymol.set_atomic_symbols_can_have_arbitrary_length(False)

    def test_mol_reader_context_accepts_mdlquiet_keyword(self):
        sdf = """ethanol
  LillyMol          2D

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.sdf")
            with open(fname, "w") as writer:
                writer.write(sdf)

            with lillymol.MolReaderContext(fname, keep_sdf_tags=True, mdlquiet=True) as reader:
                mol = reader.next()
                self.assertIsNotNone(mol)
                self.assertEqual(mol.name(), "ethanol")

            options = lillymol.ReaderOptions()
            options.keep_sdf_tags = True
            options.mdlquiet = True
            with lillymol.MolReaderContext(fname, options=options) as reader:
                mol = reader.next()
                self.assertIsNotNone(mol)
                self.assertEqual(mol.natoms(), 3)

    def test_retained_sdf_text_info_and_tags(self):
        sdf = """ethanol
  LillyMol          2D

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
>  <ID Number>
CHEMBL1

>  <Comment>
first line
second line

$$$$
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            fname = os.path.join(tmpdir, "input.sdf")
            with open(fname, "w") as writer:
                writer.write(sdf)

            with lillymol.MolReaderContext(
                    fname, lillymol.FileType.SDF, keep_sdf_tags=True) as reader:
                mol = reader.next()
                self.assertIsNotNone(mol)
                self.assertGreater(mol.number_records_text_info(), 0)
                self.assertIn(">  <ID Number>", mol.text_info())
                self.assertEqual(mol.sdf_tags()["ID_Number"], "CHEMBL1")
                self.assertEqual(mol.sdf_tags()["Comment"], "first line\nsecond line")

    def test_writer_and_context_writer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            stem = os.path.join(tmpdir, "written")
            mol = lillymol.MolFromSmiles("CCO ethanol")
            writer = lillymol.Writer()
            self.assertTrue(writer.add_output_type(lillymol.FileType.SMI))
            self.assertTrue(writer.new_stem(stem))
            self.assertTrue(writer.write(mol))
            writer.close()

            reader = lillymol.Reader()
            self.assertTrue(reader.open(stem + ".smi"))
            self.assertEqual(reader.next().name(), "ethanol")
            reader.close()

            stem = os.path.join(tmpdir, "context")
            with lillymol.MolWriterContext(stem, lillymol.FileType.SMI) as writer:
                self.assertTrue(writer.write(mol))

            reader = lillymol.Reader()
            self.assertTrue(reader.open(stem + ".smi"))
            self.assertEqual(reader.next().smiles(), "CCO")
            reader.close()

    def test_set_name(self):
        mol = lillymol.Molecule()
        self.assertTrue(mol.build_from_smiles("C methane"))
        mol.set_name("renamed")
        self.assertEqual(mol.name(), "renamed")

    def test_rdkit_style_molecule_access_aliases(self):
        mol = lillymol.MolFromSmiles("[H]OC ethanol")
        self.assertEqual(mol.GetNumAtoms(), 3)
        self.assertEqual(mol.GetNumHeavyAtoms(), 2)
        self.assertEqual(mol.GetNumBonds(), mol.nedges())
        self.assertEqual(mol.GetAtomWithIdx(1).atomic_symbol(), "O")
        self.assertEqual(mol.GetBondWithIdx(0).GetBeginAtomIdx(), 0)
        self.assertEqual([atom.atomic_symbol() for atom in mol.GetAtoms()], ["H", "O", "C"])
        self.assertEqual([(bond.a1(), bond.a2()) for bond in mol.GetBonds()], [(0, 1), (1, 2)])

        oxygen = mol.GetAtomWithIdx(1)
        mol.set_atomic_number(1, 7)
        self.assertEqual(oxygen.atomic_symbol(), "N")

    def test_molecule_common_aliases_and_counts(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertTrue(mol.ok())
        self.assertFalse(mol.empty())
        self.assertEqual(mol.GetNumAtoms(), 3)
        self.assertEqual(mol.natoms(6), 2)
        self.assertEqual(mol.natoms("O"), 1)
        with self.assertRaises(Exception):
            mol.natoms("Qq")
        self.assertEqual(mol.atomic_symbol(2), "O")
        self.assertEqual(mol.connections(1), [0, 2])
        self.assertEqual(mol.other_atom(1, 0), 0)
        self.assertEqual(mol.other_atom(1, 1), 2)

    def test_programmatic_molecule_building_and_bond_lookup(self):
        mol = lillymol.Molecule()
        self.assertTrue(mol.empty())
        c1 = mol.add_atom(6)
        c2 = mol.add_atom(6)
        o = mol.add_atom(8)
        self.assertEqual((c1, c2, o), (0, 1, 2))
        self.assertEqual(mol.natoms(), 3)
        self.assertEqual(mol.add_bond(0, 1, lillymol.BondType.SINGLE_BOND), 1)
        self.assertEqual(mol.add_bond(1, 2, lillymol.BondType.DOUBLE_BOND), 1)
        self.assertEqual(mol.nedges(), 2)
        self.assertIsNotNone(mol.bond_between_atoms(0, 1))
        self.assertIsNone(mol.bond_between_atoms(0, 2))
        self.assertEqual(mol.bond_type_between_atoms(1, 2), lillymol.BondType.DOUBLE_BOND)
        self.assertIsNone(mol.bond_type_between_atoms(0, 2))
        mol.assign_bond_numbers_to_bonds()
        self.assertTrue(mol.bond(0).bond_number_assigned())
        self.assertEqual(mol.set_bond_type_between_atoms(1, 2, lillymol.BondType.SINGLE_BOND), 1)
        self.assertEqual(mol.bond_type_between_atoms(1, 2), lillymol.BondType.SINGLE_BOND)
        self.assertEqual(mol.remove_bond_between_atoms(1, 2), 1)
        self.assertIsNone(mol.bond_between_atoms(1, 2))

        mixture = lillymol.MolFromSmiles("C methane")
        self.assertEqual(mixture.add(lillymol.MolFromSmiles("O water")), 1)
        self.assertEqual(mixture.number_fragments(), 2)

    def test_atom_removal_and_resize(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_atom(2), 1)
        self.assertEqual(mol.natoms(), 2)
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_atoms([0, 1, 0], 1), 1)
        self.assertEqual(mol.natoms(), 2)
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_atoms(lillymol.Set_of_Atoms([0, 2])), 2)
        self.assertEqual(mol.natoms(), 1)
    def test_molecule_convenience_methods(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertTrue(mol.organic_only())
        self.assertEqual(mol.non_organic_atom_count(), 0)
        self.assertTrue(mol.is_organic(0))
        self.assertTrue(mol.are_bonded(0, 1))
        self.assertFalse(mol.are_bonded(0, 2))
        self.assertEqual(mol.bonds_between(0, 2), 2)
        self.assertEqual(mol.longest_path(), 2)
        self.assertEqual(mol.most_distant_pair(), (0, 2))
        self.assertEqual(mol.atoms_on_shortest_path(0, 2), [1])
        self.assertEqual(mol.all_atoms_between(0, 2), [1])
        self.assertEqual(mol.down_the_bond(0, 1), [2])
        self.assertEqual(mol.lipinski_num_h_donors(), 1)
        self.assertEqual(mol.lipinski_num_h_acceptors(), 1)
        self.assertEqual(mol.rdkit_num_h_donors(), 1)
        self.assertEqual(mol.rdkit_num_h_acceptors(), 1)
        self.assertTrue(mol.saturated(0))
        self.assertEqual(mol.unsaturation(0), 0)

        smiles = mol.random_smiles()
        self.assertIsInstance(smiles, str)
        self.assertGreater(len(smiles), 0)
        self.assertEqual(mol.unique_kekule_smiles(), "OCC")
        self.assertEqual(mol.smiles_starting_with_atom(2), "OCC")
        order = mol.smiles_atom_order()
        self.assertEqual(sorted(order), [0, 1, 2])

        self.assertEqual(mol.renumber_atoms([2, 1, 0]), 1)
        self.assertEqual(mol.atomic_symbol(0), "O")
        with self.assertRaises(Exception):
            mol.renumber_atoms([0, 0, 1])

    def test_partial_charge_helpers(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        charges = mol.gasteiger_partial_charges()
        self.assertEqual(len(charges), mol.natoms())
        self.assertEqual(mol.partial_charge_type(), "GASTEIGER")
        self.assertAlmostEqual(mol.partial_charge(0), charges[0])
        self.assertNotEqual(sum(abs(charge) for charge in charges), 0.0)
        mol.invalidate_partial_charges()
        self.assertEqual(mol.partial_charge_type(), "")
        self.assertGreaterEqual(mol.compute_Gasteiger_partial_charges(), 0)

    def test_atoms_by_radius_single_starting_atom(self):
        mol = lillymol.MolFromSmiles("CCCCC pentane")
        shells = mol.atoms_by_radius(lillymol.Set_of_Atoms([2]), 3)
        self.assertEqual(len(shells), 4)
        self.assertCountEqual(shells[0], [2])
        self.assertCountEqual(shells[1], [1, 3])
        self.assertCountEqual(shells[2], [0, 4])
        self.assertEqual(shells[3], [])

    def test_atoms_by_radius_multiple_starting_atoms(self):
        mol = lillymol.MolFromSmiles("CCCCC pentane")
        shells = mol.atoms_by_radius(lillymol.Set_of_Atoms([0, 4]), 3)
        self.assertCountEqual(shells[0], [0, 4])
        self.assertCountEqual(shells[1], [1, 3])
        self.assertCountEqual(shells[2], [2])
        self.assertEqual(shells[3], [])

    def test_atoms_by_radius_validation(self):
        mol = lillymol.MolFromSmiles("CCC propane")
        with self.assertRaises(Exception):
            mol.atoms_by_radius(lillymol.Set_of_Atoms([3]), 1)
        with self.assertRaises(Exception):
            mol.atoms_by_radius(lillymol.Set_of_Atoms([0]), -1)

    def test_atom_map_number_helpers(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.atom_map_number(1), 0)
        mol.set_atom_map_number(1, 17)
        self.assertEqual(mol.atom_map_number(1), 17)
        self.assertEqual(mol.atom_with_atom_map_number(17), 1)
        self.assertEqual(mol.atom_with_atom_map_number(99), -1)
        mol.reset_atom_map_numbers()
        self.assertEqual(mol.atom_map_number(1), 0)

    def test_molecule_repr_debug_and_addition(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertIn("ethanol", repr(mol))
        self.assertIn("3 atoms", repr(mol))
        self.assertEqual(str(mol), "CCO ethanol")
        self.assertIn("Molecule", mol.debug_string())

        methane = lillymol.MolFromSmiles("C methane")
        combined = mol + methane
        self.assertEqual(mol.natoms(), 3)
        self.assertEqual(combined.natoms(), 4)
        mol += methane
        self.assertEqual(mol.natoms(), 4)

    def test_coordinates_object_and_transforms(self):
        coords = lillymol.Coordinates(3.0, 4.0, 0.0)
        self.assertAlmostEqual(coords.x(), 3.0)
        self.assertAlmostEqual(coords.y(), 4.0)
        self.assertAlmostEqual(coords.z(), 0.0)
        self.assertAlmostEqual(coords.norm(), 5.0)
        self.assertEqual(str(coords), "(3,4,0)")
        coords.normalise()
        self.assertAlmostEqual(coords.norm(), 1.0)
        coords.setxyz(1.0, 2.0, 3.0)
        coords.set_z(4.0)
        self.assertAlmostEqual(coords.length(), math.sqrt(21.0), delta=1.0e-6)
        self.assertAlmostEqual(coords.distance(lillymol.Coordinates(1.0, 2.0, 4.0)), 0.0)
        self.assertAlmostEqual(coords.dot_product(lillymol.Coordinates(0.0, 1.0, 0.0)), 2.0)

        mol = lillymol.MolFromSmiles("CCC propane")
        mol.set_coordinates([0.0, 0.0, 0.0,
                             1.0, 0.0, 0.0,
                             2.0, 0.0, 0.0])
        mol.translate(1.0, 2.0, 3.0)
        self.assertEqual(mol.get_coordinates(), [1.0, 2.0, 3.0,
                                                 2.0, 2.0, 3.0,
                                                 3.0, 2.0, 3.0])
        mol.translate([0, 1, 0], 1, 0.0, 1.0, 0.0)
        self.assertEqual(mol.get_coordinates()[3:6], [2.0, 3.0, 3.0])
        with self.assertRaises(Exception):
            mol.translate([1, 1], 1, 0.0, 0.0, 1.0)

        ethane = lillymol.MolFromSmiles("CC ethane")
        ethane.set_coordinates([0.0, 0.0, 0.0,
                                1.0, 0.0, 0.0])
        ethane.rotate(lillymol.Coordinates(0.0, 0.0, 1.0), math.pi / 2.0)
        self.assertAlmostEqual(ethane.x(1), 0.0, delta=1.0e-5)
        self.assertAlmostEqual(ethane.y(1), 1.0, delta=1.0e-5)

    def test_coordinates_numpy(self):
        try:
            import gc
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not available")

        mol = lillymol.MolFromSmiles("CCC propane")
        coords = np.array([0.0, 0.0, 0.0,
                           1.0, 0.0, 0.0,
                           2.0, 0.0, 0.0], dtype=np.float32)
        mol.set_coordinates_numpy(coords)

        array = mol.get_coordinates_numpy()
        self.assertEqual(array.dtype, np.dtype("float32"))
        self.assertEqual(array.shape, (9,))
        np.testing.assert_array_equal(array, coords)
        self.assertEqual(array.tolist(), mol.get_coordinates())

        updated = np.array([0.25, 0.5, 0.75,
                            1.25, 1.5, 1.75,
                            2.25, 2.5, 2.75], dtype=np.float32)
        mol.set_coordinates_numpy(updated)
        self.assertAlmostEqual(mol.x(0), 0.25)
        self.assertAlmostEqual(mol.y(1), 1.5)
        self.assertAlmostEqual(mol.z(2), 2.75)

        with self.assertRaises(Exception):
            mol.set_coordinates_numpy(updated[:-1])

        del array
        gc.collect()

    def test_dihedral_scan(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not available")

        mol = lillymol.MolFromSmiles("C{{-2,1,0}}C{{-1,0,0}}C{{0,0,0}}C{{1,1,0}} butane")
        start = mol.get_coordinates()
        conformers = mol.dihedral_scan(1, 2, 45.0)
        self.assertEqual(len(conformers), 7)
        for coords in conformers:
            self.assertEqual(coords.dtype, np.dtype("float32"))
            self.assertEqual(coords.shape, (mol.natoms() * 3,))
        self.assertEqual(mol.get_coordinates(), start)
        self.assertFalse(np.allclose(conformers[0], np.asarray(start, dtype=np.float32)))
        self.assertEqual(mol.dihedral_scan(1, 2, 180.0), [])

    def test_coordinate_and_geometry_helpers(self):
        mol = lillymol.MolFromSmiles("CCCO propanol")
        mol.setxyz(0, 0.0, 0.0, 0.0)
        mol.setxyz(1, 1.0, 0.0, 0.0)
        mol.setxyz(2, 1.0, 1.0, 0.0)
        mol.setxyz(3, 1.0, 1.0, 1.0)

        self.assertAlmostEqual(mol.x(1), 1.0)
        self.assertAlmostEqual(mol.y(2), 1.0)
        self.assertAlmostEqual(mol.z(3), 1.0)
        mol.setx(0, 0.25)
        mol.sety(0, 0.5)
        mol.setz(0, 0.75)
        self.assertEqual(mol.get_coordinates()[:3], [0.25, 0.5, 0.75])

        coords = [0.0, 0.0, 0.0,
                  1.0, 0.0, 0.0,
                  1.0, 1.0, 0.0,
                  1.0, 1.0, 1.0]
        mol.set_coordinates(coords)
        self.assertEqual(mol.get_coordinates(), coords)
        with self.assertRaises(Exception):
            mol.set_coordinates(coords[:-1])

        self.assertAlmostEqual(mol.distance_between_atoms(0, 1), 1.0)
        self.assertAlmostEqual(mol.bond_length(0, 1), 1.0)
        self.assertIsNone(mol.bond_length(0, 2))
        self.assertGreater(mol.bond_angle(1, 0, 2), 0.0)
        self.assertGreater(mol.dihedral_angle(0, 1, 2, 3), 0.0)
        self.assertNotEqual(mol.signed_dihedral_angle(0, 1, 2, 3), 0.0)
        self.assertAlmostEqual(mol.longest_intra_molecular_distance(), 3 ** 0.5)
        self.assertEqual(mol.highest_coordinate_dimensionality(), 3)
        self.assertGreater(mol.bump_check(2.0), 0)

    def test_canonical_and_symmetry_helpers(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        ranks = mol.canonical_ranks()
        self.assertEqual(len(ranks), mol.natoms())
        self.assertEqual(mol.canonical_rank(0), ranks[0])
        classes = [mol.symmetry_class(i) for i in range(mol.natoms())]
        self.assertEqual(mol.number_symmetry_classes(), len(set(classes)))

        ethane = lillymol.MolFromSmiles("CC ethane")
        equivalents = ethane.symmetry_equivalents(0)
        self.assertIn(1, equivalents)

    def test_remove_helpers(self):
        mol = lillymol.MolFromSmiles("[Na]OC sodium_ethoxide")
        self.assertEqual(mol.non_organic_atom_count(), 1)
        self.assertEqual(mol.remove_non_periodic_table_elements(), 0)

        mol = lillymol.MolFromSmiles("CC ethane")
        self.assertEqual(mol.AddHs(), 6)
        self.assertEqual(mol.remove_explicit_hydrogens(), 2)
        self.assertEqual(mol.natoms(), 2)

        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertTrue(mol.remove_bonds_to_atom(1))
        self.assertEqual(mol.nedges(), 0)

        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.remove_edge(0), 1)
        self.assertEqual(mol.nedges(), 1)
        self.assertEqual(mol.chop(1), 2)
        self.assertEqual(mol.natoms(), 2)

    def test_charges_isotopes_and_hydrogens(self):
        mol = lillymol.MolFromSmiles("C[NH3+] methylammonium")
        self.assertTrue(mol.has_formal_charges())
        self.assertEqual(mol.number_formal_charges(), 1)
        self.assertEqual(mol.net_formal_charge(), 1)
        self.assertEqual(mol.formal_charge(1), 1)
        mol.set_formal_charge(1, 0)
        self.assertEqual(mol.formal_charge(1), 0)
        self.assertEqual(mol.net_formal_charge(), 0)

        mol = lillymol.MolFromSmiles("[3CH3]CO labelled")
        self.assertEqual(mol.number_isotopic_atoms(), 1)
        self.assertEqual(mol.first_atom_with_isotope(3), 0)
        self.assertEqual(mol.transform_to_non_isotopic_form(), 1)
        self.assertEqual(mol.number_isotopic_atoms(), 0)

        mol = lillymol.MolFromSmiles("CC ethane")
        self.assertEqual(mol.implicit_hydrogens(0), 3)
        self.assertEqual(mol.explicit_hydrogens(0), 0)
        self.assertEqual(mol.hcount(0), 3)
        self.assertTrue(mol.valence_ok())
        self.assertTrue(mol.valence_ok(0))
        self.assertEqual(mol.AddHs(), 6)
        self.assertEqual(mol.natoms(), 8)
        self.assertGreater(mol.RemoveHs(), 0)
        self.assertEqual(mol.natoms(), 2)

    def test_ring_aliases_and_aromaticity(self):
        mol = lillymol.MolFromSmiles("c1ccccc1 benzene")
        self.assertTrue(mol.IsInRing(0))
        self.assertTrue(mol.in_ring_of_size(0, 6))
        self.assertTrue(mol.IsAtomInRingOfSize(0, 6))
        self.assertEqual(mol.NumAtomRings(0), 1)
        self.assertTrue(mol.is_aromatic(0))
        self.assertEqual(mol.aromatic_atom_count(), 6)
        self.assertEqual(mol.aromatic_ring_count(), 1)

    def test_chiral_centre_access_and_helpers(self):
        self.assertNotEqual(lillymol.CIP.R, lillymol.CIP.S)
        self.assertNotEqual(lillymol.CIP.Neither, lillymol.CIP.Unspecified)

        mol = lillymol.MolFromSmiles("F[C@H](Cl)Br chiral")
        self.assertEqual(mol.number_chiral_centres(), 1)
        centre = mol.chiral_centre_at_atom(1)
        self.assertIsNotNone(centre)
        self.assertEqual(centre.atom(), 1)
        self.assertTrue(centre.involves(0))
        self.assertEqual(centre.implicit_hydrogens(), 1)
        self.assertEqual(centre.lone_pairs(), 0)
        self.assertLess(lillymol.CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN, 0)
        self.assertLess(lillymol.CHIRAL_CONNECTION_IS_LONE_PAIR, 0)
        self.assertNotEqual(lillymol.CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN,
                            lillymol.CHIRAL_CONNECTION_IS_LONE_PAIR)
        self.assertEqual(centre.top_back(),
                         lillymol.CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN)
        self.assertEqual(centre.atoms(),
                         [centre.top_front(), centre.top_back(),
                          centre.left_down(), centre.right_down()])
        self.assertEqual(list(centre), centre.atoms())
        self.assertEqual([item for item in enumerate(centre)],
                         list(enumerate(centre.atoms())))
        self.assertEqual(len(centre), 4)
        self.assertEqual(centre[0], centre.top_front())
        self.assertEqual(centre[3], centre.right_down())
        with self.assertRaises(Exception):
            _ = centre[4]
        self.assertTrue(lillymol.is_chiral_implicit_hydrogen(centre.top_back()))
        self.assertFalse(lillymol.is_chiral_lone_pair(centre.top_back()))
        self.assertTrue(
            lillymol.is_chiral_lone_pair(lillymol.CHIRAL_CONNECTION_IS_LONE_PAIR))
        self.assertFalse(hasattr(centre, "atom_is_now_lone_pair"))
        self.assertFalse(hasattr(centre, "implicit_hydrogen_is_now_atom_number"))
        self.assertIn("<Chiral_Centre atom 1", repr(centre))
        centres = mol.chiral_centres()
        self.assertEqual(len(centres), 1)
        self.assertEqual(centres[0].atom(), 1)
        self.assertTrue(lillymol.is_actually_chiral(mol, 1))
        self.assertFalse(lillymol.is_actually_chiral(mol, 0))
        tag = lillymol.tetrahedral_chirality(mol, 1)
        self.assertIn(tag, [lillymol.ChiralType.CHI_TETRAHEDRAL_CW,
                            lillymol.ChiralType.CHI_TETRAHEDRAL_CCW])
        self.assertEqual(lillymol.tetrahedral_chirality(mol, 0), None)
        self.assertEqual(lillymol.tetrahedral_chirality(mol, 1, check_is_chiral=True), tag)
        self.assertEqual(mol.invert_chirality_on_atom(1), 1)
        inverted = lillymol.tetrahedral_chirality(mol, 1)
        self.assertIn(inverted, [lillymol.ChiralType.CHI_TETRAHEDRAL_CW,
                                 lillymol.ChiralType.CHI_TETRAHEDRAL_CCW])
        self.assertNotEqual(inverted, tag)
        self.assertEqual(mol.remove_chiral_centre_at_atom(1), 1)
        self.assertEqual(mol.number_chiral_centres(), 0)
        self.assertIsNone(mol.chiral_centre_at_atom(1))

    def test_remove_all_chiral_centres(self):
        mol = lillymol.MolFromSmiles("F[C@H](Cl)Br chiral")
        self.assertEqual(mol.number_chiral_centres(), 1)
        self.assertEqual(mol.remove_all_chiral_centres(), 1)
        self.assertEqual(mol.number_chiral_centres(), 0)
        self.assertIsNone(lillymol.tetrahedral_chirality(mol, 1))

    def test_atom_scalars_and_vectors(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.atomic_number(0), 6)
        self.assertEqual(mol.atomic_numbers(), [6, 6, 8])
        self.assertEqual(mol.ncon(1), 2)
        self.assertEqual(mol.nbonds(1), 2)
        self.assertEqual(mol.attached_heteroatom_count(1), 1)
        self.assertFalse(mol.is_halogen(0))
        self.assertTrue(lillymol.MolFromSmiles("Cl chloromethane").is_halogen(0))
        self.assertEqual(mol.isotopes(), [0, 0, 0])
        self.assertEqual(mol.set_isotope(1, 7), 1)
        self.assertEqual(mol.isotope(1), 7)
        self.assertEqual(mol.isotopes(), [0, 7, 0])
        self.assertEqual(mol.set_isotopes([0, 2], 5), 1)
        self.assertEqual(mol.isotopes(), [5, 7, 5])
        self.assertEqual(mol.set_isotopes(lillymol.Set_of_Atoms([1]), 9), 1)
        self.assertEqual(mol.isotopes(), [5, 9, 5])
        try:
            import numpy as np
        except ImportError:
            pass
        else:
            mol.set_isotopes(np.asarray([1, 2, 3], dtype=np.int32))
            self.assertEqual(mol.isotopes(), [1, 2, 3])
            with self.assertRaises(Exception):
                mol.set_isotopes(np.asarray([1, 2], dtype=np.int32))
        self.assertEqual(mol.remove_isotopes(), 3)
        self.assertEqual(mol.isotopes(), [0, 0, 0])

    def test_fragments_and_formula(self):
        mol = lillymol.MolFromSmiles("CC.O mixture")
        self.assertEqual(mol.number_fragments(), 2)
        self.assertEqual(mol.fragment_membership(0), 0)
        self.assertEqual(mol.get_fragment_membership(), [0, 0, 1])
        self.assertEqual(mol.atoms_in_fragment(0), 2)
        self.assertEqual(mol.atoms_in_fragment(1), 1)
        self.assertEqual(mol.atoms_in_largest_fragment(), 2)
        self.assertEqual(mol.molecular_formula(), "C2OH8")
        components = mol.create_components()
        self.assertEqual([component.natoms() for component in components], [2, 1])
        single = lillymol.MolFromSmiles("CC ethane")
        self.assertIsNone(single.create_components())

    def test_fragment_mutators(self):
        mol = lillymol.MolFromSmiles("CC.O.N mixture")
        self.assertEqual(mol.delete_fragment(1), 1)
        self.assertEqual(mol.number_fragments(), 2)
        self.assertEqual(mol.natoms(), 3)

        mol = lillymol.MolFromSmiles("CC.O.N mixture")
        self.assertEqual(mol.remove_fragment(2), 1)
        self.assertEqual(mol.natoms(), 3)

        mol = lillymol.MolFromSmiles("CC.O.N mixture")
        self.assertEqual(mol.remove_fragment_containing_atom(2), 1)
        self.assertEqual(mol.natoms(), 3)

        mol = lillymol.MolFromSmiles("CC.O mixture")
        self.assertEqual(mol.reduce_to_largest_fragment(), 1)
        self.assertEqual(mol.unique_smiles(), "CC")

        mol = lillymol.MolFromSmiles("CC.O mixture")
        self.assertEqual(mol.reduce_to_largest_fragment_carefully(), 1)
        self.assertEqual(mol.natoms(), 2)

    def test_molecular_weight(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertAlmostEqual(mol.amw(), 46.069, delta=0.001)
        self.assertAlmostEqual(lillymol.molecular_weight(mol), 46.069, delta=0.001)
        self.assertGreater(mol.exact_mass(), 0.0)

    def test_atom_access(self):
        atom = lillymol.Atom(6)
        self.assertEqual(atom.atomic_number(), 6)
        self.assertEqual(atom.atomic_symbol(), "C")
        self.assertEqual(atom.ncon(), 0)
        self.assertEqual(atom.nbonds(), 0)
        self.assertEqual(atom.formal_charge(), 0)
        self.assertTrue(atom.is_organic())
        self.assertGreater(atom.atomic_weight(), 12.0)
        self.assertEqual(atom.exact_mass(), 12.0)
        self.assertIn("<Atom C", repr(atom))

        mol = lillymol.MolFromSmiles("CCO ethanol")
        carbon = mol.atom(1)
        self.assertEqual(carbon.atomic_number(), 6)
        self.assertEqual(carbon.atomic_symbol(), "C")
        self.assertEqual(len(carbon), 2)
        self.assertEqual(carbon.connections(1), [0, 2])
        self.assertEqual([(bond.a1(), bond.a2()) for bond in carbon], [(0, 1), (1, 2)])
        self.assertEqual(carbon[0].other(1), 0)
        self.assertTrue(carbon.is_bonded_to(0))
        self.assertIn(2, carbon)
        self.assertEqual(carbon.other(1, 0), 0)
        self.assertTrue(carbon.valence_ok())
        self.assertEqual(carbon.unsaturation(), 0)
        self.assertTrue(carbon.fully_saturated())

    def test_atom_coordinate_distance_helpers(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        mol.set_coordinates([0.0, 0.0, 0.0,
                             1.0, 0.0, 0.0,
                             1.0, 1.0, 0.0])
        a0 = mol.atom(0)
        a1 = mol.atom(1)
        self.assertAlmostEqual(a0.x(), 0.0)
        self.assertAlmostEqual(a1.x(), 1.0)
        self.assertAlmostEqual(a0.distance(a1), 1.0)
        self.assertAlmostEqual(a0.distance(lillymol.Coordinates(0.0, 1.0, 0.0)), 1.0)
        self.assertAlmostEqual(a0 - a1, 1.0)
        self.assertIn("<Atom C", str(a0))

    def test_atom_view_reflects_parent_mutation(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        atom = mol.atom(1)
        self.assertEqual(atom.atomic_number(), 6)
        self.assertEqual(atom.atomic_symbol(), "C")
        self.assertEqual(mol.set_atomic_number(1, 7), 1)
        self.assertEqual(atom.atomic_number(), 7)
        self.assertEqual(atom.atomic_symbol(), "N")
        self.assertEqual(mol.atomic_number(1), 7)

    def test_bond_access(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(mol.nedges(), 2)
        bond = mol.bond(0)
        self.assertEqual(bond.a1(), 0)
        self.assertEqual(bond.a2(), 1)
        self.assertEqual(bond.other(0), 1)
        self.assertTrue(bond.involves(1))
        self.assertTrue(bond.is_single_bond())
        self.assertFalse(bond.is_double_bond())
        self.assertFalse(bond.is_triple_bond())
        self.assertFalse(bond.is_aromatic())
        self.assertIn("<Bond 0-1>", repr(bond))

        bonds = mol.bonds()
        self.assertEqual(len(bonds), 2)
        self.assertEqual([(b.a1(), b.a2()) for b in bonds], [(0, 1), (1, 2)])

    def test_bond_view_reflects_parent_mutation(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        bond = mol.bond(0)
        self.assertTrue(bond.is_single_bond())
        self.assertEqual(bond.btype(), lillymol.BondType.SINGLE_BOND)
        self.assertEqual(mol.set_bond_type_between_atoms(0, 1, lillymol.BondType.DOUBLE_BOND), 1)
        self.assertFalse(bond.is_single_bond())
        self.assertTrue(bond.is_double_bond())
        self.assertEqual(bond.btype(), lillymol.BondType.DOUBLE_BOND)

    def test_bond_rdkit_style_aliases_and_contains(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        bond = mol.bond(0)
        self.assertEqual(bond.GetBeginAtomIdx(), 0)
        self.assertEqual(bond.GetEndAtomIdx(), 1)
        self.assertEqual(bond.GetBondType(), lillymol.BondType.SINGLE_BOND)
        self.assertIn(0, bond)
        self.assertNotIn(2, bond)
        self.assertIn("<Bond 0-1>", str(bond))

    def test_bond_ring_membership(self):
        mol = lillymol.MolFromSmiles("c1ccccc1 benzene")
        bond = mol.bond(0)
        mol.ring_membership()
        self.assertEqual(bond.nrings(), 1)
        self.assertTrue(bond.IsInRing())


    def test_set_of_atoms(self):
        atoms = lillymol.Set_of_Atoms([1, 3])
        self.assertFalse(atoms.empty())
        self.assertEqual(len(atoms), 2)
        self.assertEqual(atoms.size(), 2)
        self.assertIn(1, atoms)
        self.assertNotIn(2, atoms)
        self.assertEqual(atoms[0], 1)
        self.assertEqual(list(atoms), [1, 3])
        atoms[0] = 0
        self.assertEqual(atoms.as_list(), [0, 3])
        atoms.append(4)
        atoms.extend([5, 6])
        self.assertEqual(atoms.as_list(), [0, 3, 4, 5, 6])
        values = [0] * 9
        atoms.set_vector(values, 2)
        self.assertEqual(values, [2, 0, 0, 2, 2, 2, 2, 0, 0])
        atoms.scatter(values, 4)
        self.assertEqual(values, [4, 0, 0, 4, 4, 4, 4, 0, 0])
        atoms.increment_vector(values, 1)
        self.assertEqual(values, [5, 0, 0, 5, 5, 5, 5, 0, 0])
        with self.assertRaises(Exception):
            atoms.set_vector([0, 0], 1)
        self.assertTrue(atoms.contains_both(3, 5))
        self.assertEqual(atoms.atoms_in_common(lillymol.Set_of_Atoms([5, 9])), 1)
        self.assertEqual(atoms.first_atom_in_common(lillymol.Set_of_Atoms([5, 9])), 5)
        atoms += lillymol.Set_of_Atoms([7])
        atoms += 8
        self.assertEqual(atoms.as_list()[-2:], [7, 8])
        self.assertEqual(lillymol.Set_of_Atoms([1, 2]), [1, 2])

    def test_ring_info_helpers(self):
        mol = lillymol.MolFromSmiles("c1ccccc1 benzene")
        ring_info = mol.GetRingInfo()
        self.assertEqual(ring_info.NumRings(), 1)
        self.assertEqual(ring_info.num_rings(), 1)
        self.assertEqual([sorted(ring) for ring in ring_info.AtomRings()], [[0, 1, 2, 3, 4, 5]])
        bond_rings = ring_info.BondRings()
        self.assertEqual(len(bond_rings), 1)
        self.assertEqual(len(bond_rings[0]), 6)
        self.assertEqual(ring_info.NumAtomRings(0), 1)
        self.assertEqual(ring_info.NumBondRings(0), 1)
        self.assertTrue(ring_info.IsAtomInRingOfSize(0, 6))
        self.assertTrue(ring_info.AreAtomsInSameRing(0, 3))
        self.assertTrue(ring_info.AreBondsInSameRing(0, 2))
        with self.assertRaises(Exception):
            ring_info.NumBondRings(99)

        fused = lillymol.MolFromSmiles("c1ccc2ccccc2c1 naphthalene")
        fused_info = fused.ring_info()
        self.assertEqual(fused_info.NumRings(), 2)
        self.assertEqual(len(fused_info.AtomRings()), 2)
        self.assertEqual(len(fused_info.BondRings()), 2)
        self.assertGreaterEqual(fused_info.NumAtomRings(3), 1)

    def test_ring_access(self):
        mol = lillymol.MolFromSmiles("c1ccccc1 benzene")
        self.assertEqual(mol.nrings(), 1)
        self.assertEqual(mol.get_ring_membership(), [1, 1, 1, 1, 1, 1])
        self.assertTrue(mol.is_ring_atom(0))
        self.assertEqual(mol.ring_bond_count(0), 2)
        self.assertEqual(mol.ring_bond_counts(), [2, 2, 2, 2, 2, 2])
        self.assertEqual(mol.number_ring_systems(), 1)
        self.assertEqual(mol.largest_ring_size(), 6)
        self.assertTrue(mol.in_same_ring(0, 3))
        self.assertTrue(mol.in_same_ring_system(0, 3))
        self.assertEqual(mol.fused_system_size(0), 1)

        ring = mol.ring(0)
        self.assertEqual(ring.size(), 6)
        self.assertEqual(len(ring), 6)
        self.assertEqual(ring.ring_number(), 0)
        mol.compute_aromaticity_if_needed()
        self.assertTrue(ring.is_aromatic())
        self.assertTrue(ring.contains_bond(0, 1))
        self.assertTrue(ring.contains_both(0, 3))
        self.assertFalse(ring.is_fused())
        self.assertEqual(ring.fused_ring_neighbours(), 0)
        self.assertEqual(set(ring.as_list()), set(range(6)))
        self.assertEqual(set(ring), set(range(6)))

        rings = mol.rings()
        self.assertEqual(len(rings), 1)
        self.assertEqual(rings[0].as_list(), ring.as_list())

        bridged = lillymol.MolFromSmiles("C1CC2CCC1CC2 bridged")
        self.assertEqual(bridged.nrings(), 2)
        self.assertEqual(bridged.non_sssr_rings(), 1)
        non_sssr = bridged.non_sssr_ring(0)
        self.assertEqual(non_sssr.size(), 6)
        self.assertEqual(sorted(non_sssr.as_list()), [2, 3, 4, 5, 6, 7])


    def test_molecule_sequence_and_copy(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(len(mol), 3)
        self.assertEqual(mol[2].atomic_symbol(), "O")
        self.assertEqual([atom.atomic_symbol() for atom in mol], ["C", "C", "O"])
        self.assertEqual([atom.ncon() for atom in mol], [1, 2, 1])

        lillymol.set_copy_name_in_molecule_copy_constructor(True)
        copies = [lillymol.Molecule(mol), mol.copy(), copy.copy(mol)]
        for mol_copy in copies:
            self.assertEqual(mol_copy.smiles(), mol.smiles())
            self.assertEqual(mol_copy.name(), "ethanol")
            mol_copy.set_name("copy")
            mol_copy.add_atom(6)
            self.assertEqual(mol.name(), "ethanol")
            self.assertEqual(mol.natoms(), 3)
            self.assertEqual(mol_copy.name(), "copy")
            self.assertEqual(mol_copy.natoms(), 4)

        lillymol.set_copy_name_in_molecule_copy_constructor(False)
        for unnamed_copy in (lillymol.Molecule(mol), mol.copy(), copy.copy(mol)):
            self.assertEqual(unnamed_copy.smiles(), mol.smiles())
            self.assertEqual(unnamed_copy.name(), "")
        lillymol.set_copy_name_in_molecule_copy_constructor(True)


    def test_position_3d_moves_second_fragment(self):
        mol = lillymol.MolFromSmiles("CC.CC two_ethanes")
        mol.set_coordinates([0.0, 0.0, 0.0,
                             1.5, 0.3, 0.0,
                             10.0, 0.0, 0.0,
                             11.5, -0.3, 0.0])

        self.assertEqual(lillymol.Position3D(mol, 0, 1.5, 2), 1)
        self.assertAlmostEqual(mol.distance_between_atoms(0, 2), 1.5, delta=1.0e-5)
        self.assertFalse(mol.are_bonded(0, 2))

    def test_position_3d_handles_single_atom_fragments(self):
        mol = lillymol.MolFromSmiles("C.C dimethane")
        mol.set_coordinates([0.0, 0.0, 0.0,
                             10.0, 0.0, 0.0])

        self.assertEqual(lillymol.Position3D(mol, 0, 1.5, 1), 1)
        self.assertAlmostEqual(mol.distance_between_atoms(0, 1), 1.5, delta=1.0e-5)
        self.assertFalse(mol.are_bonded(0, 1))

    def test_position_3d_rejects_atoms_in_same_fragment(self):
        mol = lillymol.MolFromSmiles("CC ethane")
        mol.set_coordinates([0.0, 0.0, 0.0,
                             1.5, 0.0, 0.0])
        self.assertEqual(lillymol.Position3D(mol, 0, 1.5, 1), 0)

    def test_molecule_contains(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertIn(6, mol)
        self.assertNotIn(7, mol)
        self.assertIn("C", mol)
        self.assertIn("O", mol)
        self.assertNotIn("N", mol)
        with self.assertRaises(Exception):
            "Qq" in mol

        carbon = lillymol.QueryFromSmarts("C")
        nitrogen = lillymol.QueryFromSmarts("N")
        self.assertIn(carbon, mol)
        self.assertNotIn(nitrogen, mol)


    def test_substructure_query_object(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        qry = lillymol.SubstructureQuery()
        self.assertTrue(qry.build_from_smarts("C"))
        self.assertEqual(qry.number_elements(), 1)
        self.assertEqual(qry.active(), 1)
        self.assertEqual(qry.substructure_search(mol), 2)
        matches = qry.substructure_search_matches(mol)
        self.assertEqual([m.as_list() for m in matches], [[0], [1]])
        self.assertEqual(qry.substructure_search_match_lists(mol), [[0], [1]])

        qry.set_max_matches_to_find(1)
        self.assertEqual(qry.substructure_search(mol), 1)

    def test_substructure_query_from_proto(self):
        mol = lillymol.MolFromSmiles("Oc1ccccc1 phenol")
        proto_string = """
query {
  smarts: "[OD1]-c:c"
  unique_embeddings_only: true
}
"""
        proto = text_format.Parse(proto_string, substructure_pb2.SubstructureQuery())
        qry = lillymol.SubstructureQuery()
        self.assertTrue(qry.construct_from_proto(proto))
        self.assertEqual(qry.substructure_search(mol), 2)
        self.assertIn(qry, mol)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".textproto", delete=False) as writer:
            writer.write(proto_string)
            fname = writer.name
        try:
            from_file = lillymol.SubstructureQuery()
            self.assertTrue(from_file.read_proto(fname))
            self.assertEqual(from_file.substructure_search(mol), 2)
        finally:
            os.remove(fname)

    def test_substructure_results_object(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        qry = lillymol.QueryFromSmarts("C")
        self.assertIsNotNone(qry)
        sresults = lillymol.SubstructureResults()
        self.assertEqual(qry.substructure_search(mol, sresults), 2)
        self.assertEqual(sresults.number_embeddings(), 2)
        self.assertEqual(len(sresults), 2)
        self.assertEqual(sresults.embeddings(), [[0], [1]])
        self.assertEqual([embedding.as_list() for embedding in sresults], [[0], [1]])
        self.assertEqual(sresults[0].as_list(), [0])
        self.assertEqual(sresults.each_embedding_set_vector(mol.natoms(), 1), [1, 1, 0])
        self.assertGreaterEqual(sresults.max_query_atoms_matched_in_search(), 1)
        sresults.reset()
        self.assertEqual(sresults.number_embeddings(), 0)

    def test_substructure_free_functions(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertTrue(lillymol.HasSubstructMatch(mol, "CO"))
        self.assertFalse(lillymol.HasSubstructMatch(mol, "N"))
        self.assertEqual(lillymol.CountSubstructMatches(mol, "C"), 2)
        self.assertEqual(lillymol.CountSubstructMatches(mol, "C", max_matches_to_find=1), 1)
        self.assertEqual(lillymol.GetSubstructMatches(mol, "CO"), [[1, 2]])
        self.assertEqual(lillymol.GetSubstructMatches(mol, "N"), [])
        self.assertIsNone(lillymol.QueryFromSmarts("["))
        with self.assertRaises(Exception):
            lillymol.HasSubstructMatch(mol, "[")


    def test_element_transformations(self):
        etrans = lillymol.ElementTransformations()
        self.assertFalse(etrans.active())
        self.assertTrue(etrans.add("I=Cl"))
        self.assertTrue(etrans.add("Br=Cl"))
        self.assertTrue(etrans.active())
        self.assertFalse(etrans.add("not a transform"))

        mol = lillymol.MolFromSmiles("ICBr halides")
        self.assertEqual(mol.natoms("I"), 1)
        self.assertEqual(mol.natoms("Br"), 1)
        self.assertEqual(etrans.process(mol), 2)
        self.assertEqual(mol.natoms("I"), 0)
        self.assertEqual(mol.natoms("Br"), 0)
        self.assertEqual(mol.natoms("Cl"), 2)

    def test_standardise(self):
        mol = lillymol.MolFromSmiles("CC(=O)[O-] acetate")
        standardise = lillymol.Standardise()
        self.assertEqual(standardise.process(mol), 0)
        self.assertEqual(mol.smiles(), "CC(=O)[O-]")

        standardise.activate_all()
        self.assertEqual(standardise.process(mol), 1)
        self.assertEqual(mol.smiles(), "CC(=O)O")

    def test_charge_assigner(self):
        charges = _charges_query_dir()
        if charges is None:
            self.skipTest("charge query directory not available")

        charge_assigner = lillymol.ChargeAssigner(charges)
        self.assertTrue(charge_assigner.active())

        mol = lillymol.MolFromSmiles("CC(=O)O acetate")
        self.assertEqual(charge_assigner.process(mol), 1)
        self.assertEqual(mol.smiles(), "CC(=O)[O-]")

        mol = lillymol.MolFromSmiles("CCN(CC)C tertiary_amine")
        self.assertEqual(charge_assigner.process(mol), 1)
        self.assertEqual(mol.smiles(), "CC[NH+](C)CC")

    def test_donor_acceptor(self):
        hbonds = _hbonds_query_dir()
        if hbonds is None:
            self.skipTest("donor/acceptor hbonds query directory not available")

        donor_acceptor = lillymol.DonorAcceptor(hbonds)
        self.assertTrue(donor_acceptor.active())

        smiles = [
            "NC1=CC=NN1 CHEMBL3217770",
            "O(C)C(=O)NN CHEMBL3183780",
            "N#CCC(=O)NN CHEMBL2106008",
        ]
        expected = [
            "[2NH2]c1[3nH][1n]cc1",
            "[1O]=C(OC)[3NH][2NH2]",
            "[1O]=C([3NH][2NH2])CC#[1N]",
        ]

        for smi, result in zip(smiles, expected):
            mol = lillymol.MolFromSmiles(smi)
            donor_acceptor.process(mol)
            self.assertEqual(mol.unique_smiles(), result)

    def test_fingerprint_default_and_tanimoto(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        bits = lillymol.linear_fingerprint(mol)
        self.assertIsInstance(bits, list)
        self.assertEqual(len(bits), 2048)
        self.assertGreater(sum(bits), 0)
        self.assertAlmostEqual(lillymol.tanimoto(bits, bits), 1.0)

        bits2 = lillymol.linear_fingerprint(mol, nbits=512, atype_specification="")
        self.assertIsNotNone(bits2)
        self.assertEqual(len(bits2), 512)
        self.assertIsNone(lillymol.linear_fingerprint(mol, nbits=512, atype_specification="BAD"))

    def test_linear_fingerprint_numpy(self):
        try:
            import gc
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not available")

        mol = lillymol.MolFromSmiles("CCO ethanol")
        bits = lillymol.linear_fingerprint(mol)
        array = lillymol.linear_fingerprint_numpy(mol)

        self.assertEqual(array.dtype, np.dtype("uint8"))
        self.assertEqual(array.shape, (2048,))
        self.assertEqual(array.tolist(), bits)
        self.assertGreater(int(array.sum()), 0)

        del array
        gc.collect()

    def test_fingerprint_creators(self):
        mol = lillymol.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine")
        ecfp = lillymol.ECFingerprintCreator(512)
        bits = ecfp.fingerprint(mol)
        self.assertIsInstance(bits, list)
        self.assertEqual(len(bits), 512)
        self.assertGreater(sum(bits), 0)

        same_molecule = lillymol.MolFromSmiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O caffeine")
        self.assertEqual(mol.unique_smiles(), same_molecule.unique_smiles())
        self.assertEqual(bits, ecfp.fingerprint(same_molecule))

        linear = lillymol.LinearFingerprintCreator(256)
        linear.set_max_length(5)
        self.assertEqual(len(linear.fingerprint(mol)), 256)

        atom_pair = lillymol.AtomPairFingerprintCreator(256)
        atom_pair.set_min_separation(1)
        atom_pair.set_max_separation(5)
        self.assertEqual(len(atom_pair.fingerprint(mol)), 256)

    def test_fingerprint_creator_numpy(self):
        try:
            import gc
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not available")

        mol = lillymol.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine")

        def assert_numpy_fingerprint_matches_list(generator, nbits):
            as_list = generator.fingerprint(mol)
            array = generator.fingerprint_numpy(mol)
            expected = np.minimum(np.asarray(as_list, dtype=np.int64), 255).astype(np.uint8)
            self.assertEqual(array.dtype, np.dtype("uint8"))
            self.assertEqual(array.shape, (nbits,))
            np.testing.assert_array_equal(array, expected)
            del array

        ecfp = lillymol.ECFingerprintCreator(512)
        assert_numpy_fingerprint_matches_list(ecfp, 512)
        ecfp.set_max_radius(1)
        assert_numpy_fingerprint_matches_list(ecfp, 512)

        linear = lillymol.LinearFingerprintCreator(256)
        linear.set_max_length(5)
        assert_numpy_fingerprint_matches_list(linear, 256)

        atom_pair = lillymol.AtomPairFingerprintCreator(256)
        atom_pair.set_min_separation(1)
        atom_pair.set_max_separation(5)
        assert_numpy_fingerprint_matches_list(atom_pair, 256)

        gc.collect()

    def test_recent_molecule_helper_methods(self):
        benzene = lillymol.MolFromSmiles("c1ccccc1 benzene")
        self.assertEqual(benzene.compute_aromaticity_if_needed(), 1)
        self.assertGreaterEqual(benzene.pi_electrons(0), 0)
        self.assertGreaterEqual(benzene.lone_pair_count(0), 0)
        self.assertIn("c", benzene.smarts_equivalent_for_atom(0))
        self.assertIn("c", benzene.smarts())
        atoms = [1] * benzene.natoms()
        self.assertTrue(benzene.find_kekule_form(atoms))
        self.assertEqual(benzene.compute_distance_matrix(), 1)
        self.assertEqual(benzene.revert_all_directional_bonds_to_non_directional(), 0)

        labelled = lillymol.MolFromSmiles("CCO ethanol")
        labelled.label_atoms_by_atom_number()
        self.assertEqual(labelled.isotopes(), [0, 1, 2])

        fused = lillymol.MolFromSmiles("C1CCC2(CC1)CCCC2 spiro")
        self.assertTrue(any(fused.is_spiro_fused(i) for i in range(fused.natoms())))
        ring_systems = fused.label_atoms_by_ring_system()
        ring_systems_spiro = fused.label_atoms_by_ring_system_including_spiro_fused()
        self.assertEqual(len(ring_systems), fused.natoms())
        self.assertEqual(len(ring_systems_spiro), fused.natoms())
        try:
            import numpy as np
        except ImportError:
            pass
        else:
            ring_systems_np = fused.label_atoms_by_ring_system_including_spiro_fused_np()
            self.assertEqual(ring_systems_np.dtype, np.dtype("int32"))
            self.assertEqual(ring_systems_np.shape, (fused.natoms(),))
            np.testing.assert_array_equal(ring_systems_np, np.asarray(ring_systems_spiro, dtype=np.int32))

        sorted_mol = lillymol.MolFromSmiles("OCN sort")
        self.assertEqual(sorted_mol.sort_atoms([2, 1, 0]), 1)
        self.assertEqual(sorted_mol.atomic_numbers(), [8, 6, 7])
        with self.assertRaises(Exception):
            sorted_mol.sort_atoms([1, 2])

        moved = lillymol.MolFromSmiles("OCC move")
        self.assertEqual(moved.move_to_end_of_connection_table(8), 1)
        self.assertEqual(moved.atomic_numbers()[-1], 8)

        scaffold_source = lillymol.MolFromSmiles("CCc1ccccc1C(=O)O scaffold")
        scaffold = scaffold_source.scaffold()
        self.assertLess(scaffold.natoms(), scaffold_source.natoms())
        self.assertGreater(scaffold_source.to_scaffold(), 0)
        self.assertEqual(scaffold_source.unique_smiles(), scaffold.unique_smiles())

        graph = lillymol.MolFromSmiles("CCO graph")
        self.assertEqual(graph.change_to_graph_form(), 1)
        graph2 = lillymol.MolFromSmiles("CCO graph")
        mol2graph = lillymol.Mol2Graph()
        mol2graph.turn_on_most_useful_options()
        self.assertEqual(graph2.to_graph(mol2graph), 1)

    def test_element_and_hybridization_helpers(self):
        self.assertEqual(lillymol.count_atoms_in_smiles("CCO"), 3)
        self.assertEqual(lillymol.count_atoms_in_smiles("c1ccccc1"), 6)

        mol = lillymol.MolFromSmiles("CC#N acetonitrile")
        self.assertEqual(lillymol.hybridization(mol, 0), lillymol.Hybridization.SP3)
        self.assertEqual(lillymol.hybridization(mol, 1), lillymol.Hybridization.SP)
        self.assertEqual(mol.hybridization(1), lillymol.Hybridization.SP)
        self.assertEqual(lillymol.hybridization_name(lillymol.Hybridization.SP3), "SP3")
        with self.assertRaises(Exception):
            lillymol.hybridization(mol, 99)
        with self.assertRaises(Exception):
            mol.hybridization(99)

        lillymol.set_auto_create_new_elements(0)
        lillymol.set_atomic_symbols_can_have_arbitrary_length(0)
        lillymol.set_display_strange_chemistry_messages(1)
        lillymol.set_display_smiles_interpretation_error_messages(1)
        self.assertIn(lillymol.interpret_D_as_deuterium(), [0, 1])
        self.assertIn(lillymol.interpret_T_as_deuterium(), [0, 1])

    def test_mformula_build_counts_and_fingerprint(self):
        formula = lillymol.MFormula()
        self.assertFalse(formula.initialised())
        self.assertEqual(formula.build_from_smiles("CCO"), 3)
        self.assertTrue(formula.initialised())
        self.assertEqual(formula.natoms(), 3)
        self.assertEqual(formula.carbon(), 2)
        self.assertEqual(formula.oxygen(), 1)
        self.assertEqual(formula.nitrogen(), 0)
        self.assertEqual(len(formula.fixed_counted_fingerprint()), 18)

    def test_mformula_build_from_molecule_and_subset(self):
        ethanol = lillymol.MolFromSmiles("CCO ethanol")
        ethoxy = lillymol.MolFromSmiles("CCOC ethoxy")
        formula = lillymol.MFormula()
        larger = lillymol.MFormula()

        self.assertEqual(formula.build(ethanol), 1)
        self.assertEqual(larger.build(ethoxy), 1)
        self.assertEqual(formula.carbon(), 2)
        self.assertEqual(formula.oxygen(), 1)
        self.assertEqual(larger.carbon(), 3)
        self.assertEqual(larger.oxygen(), 1)
        self.assertTrue(formula.is_subset(formula))
        self.assertTrue(formula.is_element_count_subset(larger))
        self.assertGreater(larger.diff(formula), 0)

    def test_mformula_build_from_selected_atoms(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        formula = lillymol.MFormula()
        self.assertEqual(formula.build(mol, lillymol.Set_of_Atoms([1, 2])), 1)
        self.assertEqual(formula.carbon(), 1)
        self.assertEqual(formula.oxygen(), 1)

    def test_atom_typing_specification(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        atom_typing = lillymol.AtomTypingSpecification("UST:Y")
        self.assertTrue(atom_typing.active())

        types = atom_typing.assign_atom_types(mol)
        self.assertEqual(len(types), mol.natoms())
        self.assertEqual(types, lillymol.assign_atom_types(mol, "UST:Y"))
        self.assertEqual(types[0], types[1])
        self.assertNotEqual(types[1], types[2])
        self.assertTrue(all(isinstance(value, int) for value in types))
        with self.assertRaises(Exception):
            atom_typing.string_representation()
        self.assertTrue(atom_typing.append_to_tag("FP").startswith("FP"))

        atomic_number_typing = lillymol.AtomTypingSpecification("z")
        self.assertEqual(atomic_number_typing.string_representation(), "z")

        empty = lillymol.AtomTypingSpecification()
        self.assertFalse(empty.active())
        self.assertTrue(empty.build("UST:AY"))
        self.assertTrue(empty.active())
        with self.assertRaises(Exception):
            lillymol.AtomTypingSpecification("BAD")

    def test_alogp_class(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        calc = lillymol.ALogP()
        self.assertAlmostEqual(calc.logp(mol), lillymol.alogp(mol), delta=1.0e-6)
        calc.set_rdkit_phoshoric_acid_hydrogen(True)
        calc.set_use_alcohol_for_acid(True)
        self.assertIsNotNone(calc.logp(mol))

    def test_tpsa_class(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        tpsa = lillymol.TPSA()
        value = tpsa.compute(mol)
        self.assertIsNotNone(value)
        self.assertIsInstance(value, float)
        self.assertGreater(value, 0.0)
        self.assertAlmostEqual(value, lillymol.tpsa(mol), delta=1.0e-6)

    def test_tpsa_empty_molecule_returns_none(self):
        self.assertIsNone(lillymol.TPSA().compute(lillymol.Molecule()))

    def test_tpsa_options_and_rdkit_compatibility(self):
        tpsa = lillymol.TPSA()
        tpsa.set_display_psa_unclassified_atom_messages(0)
        tpsa.set_return_zero_for_unclassified_atoms(0)
        tpsa.set_non_zero_contribution_for_SD2(1)
        tpsa.set_zero_for_all_sulphur_atoms(0)
        tpsa.set_zero_for_all_phosphorus_atoms(0)
        tpsa.set_convert_to_charge_separated(0)
        self.assertEqual(tpsa.display_psa_unclassified_atom_messages(), 0)
        self.assertEqual(tpsa.zero_for_all_sulphur_atoms(), 0)

        tpsa.set_rdkit_compatibility()
        self.assertEqual(tpsa.zero_for_all_sulphur_atoms(), 1)
        self.assertEqual(tpsa.zero_for_all_phosphorus_atoms(), 1)
        self.assertEqual(tpsa.convert_to_charge_separated(), 1)
        self.assertIsNotNone(tpsa.compute(lillymol.MolFromSmiles("CCO ethanol")))

    def _medchemwizard(self):
        reactions = _medchemwizard_reactions_file()
        if reactions is None:
            self.skipTest("Cannot locate MedchemWizard REACTIONS data")

        wizard = lillymol.MedchemWizard()
        wizard.read_reactions(reactions)
        wizard.set_max_atoms(3)
        wizard.set_append_names(True)
        wizard.set_name_separator(" ")
        return wizard

    def test_medchemwizard_read_single_reaction_with_relative_query(self):
        with tempfile.TemporaryDirectory() as directory:
            query = os.path.join(directory, "query.textproto")
            with open(query, "w") as output:
                output.write('query { smarts: "[C]-[OD1]" }')
            reaction = os.path.join(directory, "reaction.textproto")
            with open(reaction, "w") as output:
                output.write('scaffold { query_file: "query.textproto" '
                             'change_element { atom: 1 element: "N" } }')
            wizard = lillymol.MedchemWizard()
            wizard.read_reaction(reaction)
            mol = lillymol.MolFromSmiles("CO methanol")
            products = wizard.process(mol)
            self.assertEqual(len(products), 1)
            self.assertEqual(products[0].unique_smiles(),
                             lillymol.MolFromSmiles("CN").unique_smiles())
            with self.assertRaises(RuntimeError):
                wizard.read_reaction(os.path.join(directory, "missing.rxn"))
            self.assertEqual(wizard.number_reactions(), 1)
            wizard.read_reaction(reaction)
            self.assertEqual(wizard.number_reactions(), 2)
            self.assertEqual(len(wizard.process(mol)), 2)

    def test_medchemwizard_initialise_from_environment(self):
        if "LILLYMOL_HOME" not in os.environ:
            self.skipTest("LILLYMOL_HOME not set")

        wizard = lillymol.MedchemWizard()
        wizard.initialise_from_environment()
        self.assertGreater(wizard.number_reactions(), 0)

    def test_medchemwizard_products_do_not_change_input(self):
        wizard = self._medchemwizard()
        mol = _methane()
        initial_smiles = mol.smiles()

        products = wizard.process(mol)

        self.assertGreater(len(products), 0)
        self.assertEqual(mol.smiles(), initial_smiles)
        self.assertTrue(all(isinstance(product, lillymol.Molecule)
                            for product in products))
        self.assertTrue(all(product.natoms() <= 3 for product in products))

        stats = wizard.stats()
        self.assertEqual(stats["molecules_read"], 1)
        self.assertGreaterEqual(stats["molecules_produced"], len(products))

    def test_medchemwizard_protected_atoms_suppress_products(self):
        wizard = self._medchemwizard()
        wizard.add_do_not_change_smarts("[*]")

        products = wizard.process(_methane())

        self.assertEqual(len(products), 0)
        self.assertGreater(
            wizard.stats()["embeddings_rejected_for_changing_protected_atoms"], 0)

    def test_medchemwizard_protected_query_no_match_can_be_ignored(self):
        wizard = self._medchemwizard()
        wizard.add_do_not_change_smarts("[Cl]")

        with self.assertRaises(RuntimeError):
            wizard.process(_methane())

        wizard = self._medchemwizard()
        wizard.add_do_not_change_smarts("[Cl]")
        wizard.set_ignore_do_not_change_queries_not_matching(True)
        self.assertGreater(len(wizard.process(_methane())), 0)

    def test_dicer_default_and_break_cc(self):
        dicer = lillymol.Dicer()
        mol = lillymol.MolFromSmiles("CCCC butane")
        self.assertEqual(dicer.dice(mol), {})

        dicer.set_break_cc_bonds(True)
        self.assertEqual(dicer.dice(mol), {"CC": 2, "CCC": 2, "C": 2})

    def test_dicer_join_point_labels_and_size_limits(self):
        dicer = lillymol.Dicer()
        dicer.set_label_join_points(8)
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertEqual(dicer.dice(mol), {"[8CH3]C": 1, "[8OH2]": 1})

        dicer = lillymol.Dicer()
        dicer.set_label_join_points(8)
        dicer.set_min_fragment_size(2)
        self.assertEqual(dicer.dice(mol), {"[8CH3]C": 1})

        dicer = lillymol.Dicer()
        dicer.set_label_join_points(8)
        dicer.set_max_fragment_size(1)
        self.assertEqual(dicer.dice(mol), {"[8OH2]": 1})

    def test_dicer_increment_join_points_and_global_counts(self):
        dicer = lillymol.Dicer()
        dicer.set_increment_isotope_for_join_points(100)
        mol = lillymol.MolFromSmiles("CC[23CH2]NCC")
        self.assertEqual(dicer.dice(mol), {
            "[100CH3]C": 1,
            "[100NH2]CC": 1,
            "[100NH2][100CH2]CC": 1,
            "[123CH3]CC": 1,
        })

        dicer = lillymol.Dicer()
        dicer.set_accumulate_global_fragment_count(True)
        dicer.set_break_cc_bonds(True)
        self.assertEqual(dicer.dice(lillymol.MolFromSmiles("CC ethane")), {"C": 2})
        self.assertEqual(dicer.get_global_fragment_count(), {"C": 1})
        self.assertEqual(dicer.dice(lillymol.MolFromSmiles("CC ethane")), {"C": 2})
        self.assertEqual(dicer.get_global_fragment_count(), {"C": 2})

    def test_dicer_break_bond_smarts(self):
        dicer = lillymol.Dicer()
        dicer.set_label_join_points(8)
        dicer.set_break_cc_bonds(True)
        dicer.set_max_bonds_to_break(3)
        self.assertTrue(dicer.add_bond_break_smarts("C-F"))
        self.assertFalse(dicer.add_bond_break_smarts("["))
        mol = lillymol.MolFromSmiles("CC(F)(F)F trifluoroethane")
        self.assertEqual(dicer.dice(mol), {
            "[8FH]": 15,
            "[8CH3]C": 6,
            "F[8CH2]C": 6,
            "C[8CH](F)F": 3,
        })

    def test_ring_replacement_no_match(self):
        replacement_file = _ring_replacement_file()
        if replacement_file is None:
            self.skipTest("ring replacement test data not available")

        replacement = lillymol.RingReplacement()
        self.assertTrue(replacement.set_ring_atom_smarts("[#7]"))
        self.assertGreater(replacement.read_replacement_rings(replacement_file), 0)
        self.assertGreater(replacement.number_replacement_rings(), 0)

        mol = lillymol.MolFromSmiles("Oc1ccc(OC)cc1")
        self.assertEqual(replacement.process(mol), [])

    def test_ring_replacement_products(self):
        replacement_file = _ring_replacement_file()
        if replacement_file is None:
            self.skipTest("ring replacement test data not available")

        replacement = lillymol.RingReplacement()
        self.assertTrue(replacement.set_ring_atom_smarts("[1c]"))
        self.assertGreater(replacement.read_replacement_rings(replacement_file), 0)

        mol = lillymol.MolFromSmiles("O[1c]1cc[1c](OC)cc1")
        products = replacement.process(mol)
        self.assertEqual(len(products), 12)
        self.assertTrue(all(isinstance(product, lillymol.Molecule)
                            for product in products))
        self.assertTrue(all(product.number_isotopic_atoms() == 2
                            for product in products))

    def test_ring_replacement_unique_only(self):
        replacement_file = _ring_replacement_file()
        if replacement_file is None:
            self.skipTest("ring replacement test data not available")

        replacement = lillymol.RingReplacement()
        self.assertTrue(replacement.set_ring_atom_smarts("c"))
        replacement.set_unique_molecules_only(True)
        replacement.set_min_support_requirement(100)
        self.assertGreater(replacement.read_replacement_rings(replacement_file), 0)

        mol = lillymol.MolFromSmiles("Oc1ccc(OC)cc1 start")
        self.assertEqual(len(replacement.process(mol)), 7)
        self.assertEqual(replacement.process(mol), [])
        replacement.clear_unique_molecule_cache()
        self.assertEqual(len(replacement.process(mol)), 7)

    def test_ring_replacement_remove_isotopes(self):
        replacement_file = _ring_replacement_file()
        if replacement_file is None:
            self.skipTest("ring replacement test data not available")

        replacement = lillymol.RingReplacement()
        self.assertTrue(replacement.set_ring_atom_smarts("[1c]"))
        replacement.set_remove_isotopes(True)
        self.assertGreater(replacement.read_replacement_rings(replacement_file), 0)

        mol = lillymol.MolFromSmiles("O[1c]1cc[1c](OC)cc1")
        products = replacement.process(mol)
        self.assertEqual(len(products), 12)
        self.assertTrue(all(product.number_isotopic_atoms() == 0
                            for product in products))

    def test_ring_replacement_preserve_connectivity(self):
        replacement_file = os.path.join(
            os.path.dirname(__file__), "testdata",
            "ring_replacement_connectivity.textproto")

        replacement = lillymol.RingReplacement()
        self.assertEqual(replacement.read_replacement_rings(replacement_file), 1)
        with self.assertRaises(RuntimeError):
            replacement.set_preserve_connectivity(False)

        replacement = lillymol.RingReplacement()
        replacement.set_preserve_connectivity(False)
        self.assertEqual(replacement.read_replacement_rings(replacement_file), 2)

    def test_reaction_rdkit_cookbook_smirks(self):
        core = lillymol.MolFromSmiles("*c1c(C)cccc1O")
        sidechain = lillymol.MolFromSmiles("CN*")

        lillymol.set_smirks_lost_atom_means_remove_frgment(1)

        reaction = lillymol.Reaction()
        self.assertTrue(reaction.construct_from_smirks(
            "[c:1][#0:3].[#0:4][*:2]>>[*:1]-[*:2]"))
        smc = lillymol.SidechainMatchConditions()
        self.assertTrue(reaction.add_sidechain_reagent(0, sidechain, smc))
        self.assertEqual(reaction.number_sidechains(), 1)
        self.assertEqual(reaction.number_sidechains_with_reagents(), 1)
        self.assertEqual(reaction.sidechain_name(0, 0), sidechain.name())

        products = reaction.perform_reaction(core, sidechain)
        self.assertIsNotNone(products)
        self.assertEqual(len(products), 1)
        self.assertEqual(products[0].unique_smiles(), "Oc1c(NC)c(C)ccc1")

    def test_reaction_textproto_multiple_reagents(self):
        reagents = [
            lillymol.MolFromSmiles("O-C(=O)c1ccc(Cl)cc1 scaffold"),
            lillymol.MolFromSmiles("Nc1ccc(S)cc1 R1"),
            lillymol.MolFromSmiles("C R2"),
        ]
        reaction_textproto = """scaffold {
  id: 0
  smarts: "[OD1]-C=O.[Cl:3]"
  remove_atom: 0
  remove_atom: 3
}
sidechain {
  id: 1
  smarts: "[N]-[c:1].[c:2]-[S:3]"
  remove_atom: 3
  join {
    a1: 1
    a2: 0
  }
}
sidechain {
  id: 2
  smarts: "C"
  join {
    c1 {
      component: 1
      atom: 2
    }
    a2: 0
  }
}
"""
        reaction = lillymol.Reaction()
        self.assertTrue(reaction.construct_from_textproto(reaction_textproto))
        self.assertEqual(reaction.number_sidechains(), 2)

        product = reaction.perform_reaction(reagents)
        self.assertIsNotNone(product)
        self.assertEqual(product.unique_smiles(), "O=C(Nc1ccc(C)cc1)c1ccccc1")

    def test_reaction_iterator_and_reagent_names(self):
        core = lillymol.MolFromSmiles("*c1ccccc1 scaffold")
        sidechain = lillymol.MolFromSmiles("CN* methylamine")

        reaction = lillymol.Reaction()
        self.assertTrue(reaction.construct_from_smirks(
            "[c:1][#0:3].[#0:4][*:2]>>[*:1]-[*:2]"))
        smc = lillymol.SidechainMatchConditions()
        self.assertTrue(reaction.add_sidechain_reagent(0, sidechain, smc))

        iterator = lillymol.ReactionIterator(reaction)
        self.assertTrue(iterator.active())
        self.assertEqual(iterator.reagent(0), 0)
        self.assertEqual(reaction.reagent_names(iterator), ["methylamine"])

        matches = reaction.substructure_search_matches(core)
        self.assertIsNotNone(matches)
        product = reaction.perform_reaction(core, lillymol.Set_of_Atoms(matches[0]), iterator)
        self.assertIsNotNone(product)
        self.assertEqual(product.unique_smiles(), "CNc1ccccc1")

    def test_gfp_list_read_and_distance(self):
        fname = _write_text_file(CARBON_GFP, suffix=".gfp")
        try:
            gfp = lillymol.GFPList.from_file(fname)
            self.assertEqual(len(gfp), 2)
            self.assertEqual(gfp.size(), 2)
            self.assertEqual(gfp.tags(), ["FCTS<"])
            self.assertEqual(gfp.smiles(0), "C")
            self.assertEqual(gfp.id(1), "Ethane")
            self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
            self.assertAlmostEqual(gfp.distance(1, 0), 0.5, places=6)
        finally:
            os.remove(fname)

    def test_gfp_nearest_neighbours(self):
        fname = _write_text_file(CARBON_GFP, suffix=".gfp")
        try:
            gfp = lillymol.GFPList.from_file(fname)
            hits = gfp.nearest_neighbours(0, 1)
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0].index, 1)
            self.assertAlmostEqual(hits[0].distance, 0.5, places=6)
            self.assertIn("GFPNearestNeighbour", repr(hits[0]))

            close = gfp.nearest_neighbours_within_distance(0, 0.5)
            self.assertEqual(len(close), 1)
            self.assertEqual(close[0].index, 1)
            self.assertEqual(gfp.nearest_neighbours_within_distance(0, 0.49), [])
        finally:
            os.remove(fname)

    def test_gfp_standard_golden_distances(self):
        gfp = lillymol.GFPList.from_file(_testdata_file("rand10.standard.gfp"))
        self.assertGreater(len(gfp), 3)
        self.assertEqual(gfp.id(0), "CHEMBL3460651")
        self.assertEqual(gfp.id(1), "CHEMBL3460651.a")
        self.assertEqual(gfp.id(3), "CHEMBL1417367")

        self.assertAlmostEqual(gfp.distance(0, 1), 0.0421, delta=0.0001)
        self.assertAlmostEqual(gfp.distance(1, 0), 0.0421, delta=0.0001)
        self.assertAlmostEqual(gfp.distance(3, 0), 0.499, delta=0.001)
        self.assertAlmostEqual(gfp.distance(0, 3), 0.499, delta=0.001)

    def test_gfp_generator_specs_match_standard_context(self):
        standard = lillymol.GFPContext.standard()
        from_specs = lillymol.GFPContext.from_specs([
            lillymol.GFP.iw(),
            lillymol.GFP.maccs(),
            lillymol.GFP.mpr(),
        ])

        self.assertEqual(from_specs.tags(), standard.tags())
        mol = lillymol.MolFromSmiles("CCO ethanol")
        fp_standard = standard.fingerprint(mol)
        fp_from_specs = from_specs.fingerprint(mol)
        self.assertEqual(fp_standard.context_hash(), fp_from_specs.context_hash())
        self.assertAlmostEqual(standard.distance(fp_standard, fp_from_specs),
                               0.0, places=6)

    def test_gfp_standard_list_add_and_query_fingerprint(self):
        gfp = lillymol.GFPList.standard()
        self.assertEqual(gfp.tags(), ["FPIW<", "FPMK<", "FPMK2<", "MPR<"])

        for smiles in ["CC ethane", "CCC propane", "CCCC butane"]:
            gfp.add(lillymol.MolFromSmiles(smiles))

        self.assertEqual(len(gfp), 3)
        self.assertEqual(gfp.id(0), "ethane")
        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)

        query = lillymol.MolFromSmiles("CCC query")
        fp = lillymol.GFPContext.standard().fingerprint(query)
        self.assertAlmostEqual(gfp.distance(fp, 1), 0.0, places=6)

        hits = gfp.nearest_neighbours(fp, 2)
        self.assertEqual(hits[0].index, 1)
        self.assertAlmostEqual(hits[0].distance, 0.0, places=6)

    def test_gfp_programmatic_standard_list_matches_stored_file_distances(self):
        stored = lillymol.GFPList.from_file(_testdata_file("rand10.standard.gfp"))
        molecules = [lillymol.MolFromSmiles(stored.smiles(i))
                     for i in range(len(stored))]
        generated = lillymol.GFPList.standard_from_molecules(molecules)

        self.assertEqual(len(generated), len(stored))
        self.assertFalse(generated.metadata_stored())
        with self.assertRaises(RuntimeError):
            generated.id(0)
        with self.assertRaises(RuntimeError):
            generated.smiles(0)

        for i in range(len(stored)):
            for j in range(len(stored)):
                self.assertAlmostEqual(generated.distance(i, j),
                                       stored.distance(i, j),
                                       places=6,
                                       msg=f"{i},{j}")

    def test_gfp_standard_from_molecules_can_store_metadata(self):
        molecules = [lillymol.MolFromSmiles("CC ethane"),
                     lillymol.MolFromSmiles("CCC propane"),
                     lillymol.MolFromSmiles("CCCC butane")]
        gfp = lillymol.GFPList.standard_from_molecules(molecules, store_metadata=True)

        self.assertEqual(len(gfp), 3)
        self.assertTrue(gfp.metadata_stored())
        self.assertEqual(gfp.id(0), "ethane")
        self.assertEqual(gfp.smiles(1), "CCC")
        self.assertAlmostEqual(gfp.distance(1, 1), 0.0, places=6)

    def test_gfp_add_molecules_can_store_metadata(self):
        molecules = [lillymol.MolFromSmiles("CC ethane"),
                     lillymol.MolFromSmiles("CCC propane")]
        gfp = lillymol.GFPList.standard()
        gfp.add_molecules(molecules, store_metadata=True)

        self.assertEqual(len(gfp), 2)
        self.assertEqual(gfp.id(1), "propane")
        self.assertEqual(gfp.smiles(0), "CC")

    def test_gfp_add_smiles_matches_add_molecules_without_metadata(self):
        smiles = ["CC", "CCC", "CCCC"]
        from_smiles = lillymol.GFPList.standard()
        from_smiles.add_smiles(smiles)

        molecules = [lillymol.MolFromSmiles(smi) for smi in smiles]
        from_molecules = lillymol.GFPList.standard_from_molecules(molecules)

        self.assertEqual(len(from_smiles), len(smiles))
        self.assertFalse(from_smiles.metadata_stored())
        with self.assertRaises(RuntimeError):
            from_smiles.id(0)

        for i in range(len(smiles)):
            for j in range(len(smiles)):
                self.assertAlmostEqual(from_smiles.distance(i, j),
                                       from_molecules.distance(i, j),
                                       places=6)

    def test_gfp_add_smiles_rejects_bad_smiles(self):
        gfp = lillymol.GFPList.standard()
        with self.assertRaisesRegex(RuntimeError, "index 1"):
            gfp.add_smiles(["CC", "["])

    def test_gfp_add_smiles_rejects_existing_metadata_mode(self):
        gfp = lillymol.GFPList.standard()
        gfp.add(lillymol.MolFromSmiles("CC ethane"))
        with self.assertRaises(RuntimeError):
            gfp.add_smiles(["CCC"])

    def test_gfp_generator_specs_and_errors(self):
        self.assertEqual(lillymol.GFP.maccs(level2=False).components(), ["FPMK<"])
        self.assertEqual(lillymol.GFP.formula().components(), ["FCFML<"])
        self.assertEqual(lillymol.GFP.atom_pair(min_separation=0, max_separation=2,
                                                   include_out_of_range=True).components(),
                         ["NCAPT0M2USTY<"])
        self.assertEqual(lillymol.GFP.ec(radius=3, atom_type="UST:AY").components(),
                         ["NCEC3USTAY<"])
        self.assertEqual(lillymol.GFP.spinach(label_join_points=True).components(),
                         ["FPSPINI<"])
        self.assertEqual(lillymol.GFP.scaffold(label_join_points=True).components(),
                         ["FPSCAFI<"])
        substructure = lillymol.GFP.substructure("c1ccccc1", radius=1)
        self.assertTrue(substructure.components()[0].startswith("FPSUB1USTARY"))

        with self.assertRaises(ValueError):
            lillymol.GFP.alogp(replicates=0)
        with self.assertRaises(ValueError):
            lillymol.GFP.atom_pair(min_separation=4, max_separation=3)
        with self.assertRaises(ValueError):
            lillymol.GFP.substructure("C", no_match="skip")
        with self.assertRaises(RuntimeError):
            lillymol.GFPContext.from_specs([lillymol.GFP.iw(), lillymol.GFP.iw()])

    def _run_truncated_distance_matrix_storage_case(self, storage):
        fname = _write_tfdatarecord([
            _nearneighbours("A", [("B", 0.25), ("C", 0.50)]),
            _nearneighbours("B", [("A", 0.25)]),
            _nearneighbours("C", [("A", 0.50)]),
            _nearneighbours("D", []),
        ])
        try:
            dm = lillymol.TruncatedDistanceMatrix(fname, storage=storage)
            self.assertEqual(len(dm), 4)
            self.assertEqual(dm.size(), 4)
            self.assertEqual(dm.number_distances(), 2)
            self.assertEqual(dm.index("A"), 0)
            self.assertIsNone(dm.index("missing"))
            self.assertEqual(dm.name(2), "C")

            self.assertAlmostEqual(dm.distance(0, 1), 0.25, delta=1.0 / 255.0)
            self.assertAlmostEqual(dm.distance(1, 0), 0.25, delta=1.0 / 255.0)
            self.assertAlmostEqual(dm.distance(0, 2), 0.50, delta=1.0 / 255.0)
            self.assertIsNone(dm.distance(1, 2))
            self.assertAlmostEqual(dm.distance_or_default(1, 2),
                                   dm.max_stored_distance(), places=6)
            self.assertAlmostEqual(dm.distance_or_default(3, 3), 0.0, places=6)

            with self.assertRaises(ValueError):
                dm.set_default_distance(0.25)
            dm.set_default_distance(1.0)
            self.assertAlmostEqual(dm.distance_or_default(1, 2), 1.0, places=6)
            self.assertEqual(dm.distances_or_default([0, 1, 3], [1, 2, 3]),
                             [dm.distance(0, 1), 1.0, 0.0])
            with self.assertRaises(Exception):
                dm.name(99)
            with self.assertRaises(Exception):
                dm.distances_or_default([0], [1, 2])
        finally:
            os.remove(fname)

    def test_truncated_distance_matrix_row_sparse(self):
        self._run_truncated_distance_matrix_storage_case(
            lillymol.TruncatedDistanceMatrixStorage.ROW_SPARSE)

    def test_truncated_distance_matrix_row_hash(self):
        self._run_truncated_distance_matrix_storage_case(
            lillymol.TruncatedDistanceMatrixStorage.ROW_HASH)

    def test_truncated_distance_matrix_indexed_proto(self):
        fname = _write_tfdatarecord([
            _nearneighbours_indices("A", [(1, 0.25), (2, 0.50)]),
            _nearneighbours_indices("B", [(0, 0.25)]),
            _nearneighbours_indices("C", [(0, 0.50)]),
            _nearneighbours_indices("D", []),
        ])
        try:
            dm = lillymol.TruncatedDistanceMatrix(
                fname,
                storage=lillymol.TruncatedDistanceMatrixStorage.ROW_SPARSE,
                proto_type=lillymol.TruncatedDistanceMatrixProto.NEARNEIGHBOURS_INDICES)
            self.assertEqual(dm.size(), 4)
            self.assertEqual(dm.number_distances(), 2)
            self.assertEqual(dm.name(1), "B")
            self.assertAlmostEqual(dm.distance(2, 0), 0.50, delta=1.0 / 255.0)
        finally:
            os.remove(fname)

    def test_truncated_distance_matrix_one_byte_duplicate(self):
        fname = _write_tfdatarecord([
            _nearneighbours("A", [("B", 10.0 / 255.0)]),
            _nearneighbours("B", [("A", 11.0 / 255.0)]),
        ])
        try:
            dm = lillymol.TruncatedDistanceMatrix(fname)
            self.assertEqual(dm.number_distances(), 1)
            self.assertEqual(dm.duplicate_distances_differing_by_one(), 1)
            self.assertAlmostEqual(dm.distance(0, 1), 10.0 / 255.0, places=6)
        finally:
            os.remove(fname)

    def test_truncated_distance_matrix_rejects_conflicting_duplicate(self):
        fname = _write_tfdatarecord([
            _nearneighbours("A", [("B", 0.25)]),
            _nearneighbours("B", [("A", 0.75)]),
        ])
        try:
            with self.assertRaises(RuntimeError):
                lillymol.TruncatedDistanceMatrix(fname)
        finally:
            os.remove(fname)

    def test_unique_molecules(self):
        unique = lillymol.UniqueMolecules()
        self.assertTrue(unique.is_unique(lillymol.MolFromSmiles("C methane")))
        self.assertFalse(unique.is_unique(lillymol.MolFromSmiles("C methane_again")))

        chiral = lillymol.UniqueMolecules()
        self.assertTrue(chiral.is_unique(lillymol.MolFromSmiles("C(O)[C@@H](N)C")))
        self.assertTrue(chiral.is_unique(lillymol.MolFromSmiles("C(O)[C@H](N)C")))

        achiral = lillymol.UniqueMolecules()
        achiral.set_include_chiral_info(False)
        self.assertTrue(achiral.is_unique(lillymol.MolFromSmiles("C(O)[C@@H](N)C")))
        self.assertFalse(achiral.is_unique(lillymol.MolFromSmiles("C(O)[C@H](N)C")))

        fragments = lillymol.UniqueMolecules()
        fragments.set_strip_to_largest_fragment(True)
        self.assertTrue(fragments.is_unique(lillymol.MolFromSmiles("CC.C")))
        self.assertFalse(fragments.is_unique(lillymol.MolFromSmiles("CC.O")))

        isotopes = lillymol.UniqueMolecules()
        isotopes.set_consider_isotopes(False)
        mol = lillymol.MolFromSmiles("C methane")
        self.assertTrue(isotopes.is_unique(mol))
        mol.set_isotope(0, 1)
        self.assertFalse(isotopes.is_unique(mol))

        graph = unique.graph_specifications()
        self.assertFalse(graph.active())
        graph.turn_on_most_useful_options()
        self.assertTrue(graph.active())
        graph.set_active(False)
        self.assertFalse(graph.active())
        self.assertEqual(unique.report(), 1)

    def test_qed(self):
        query_dir = _qed_query_dir()
        if query_dir is None:
            self.skipTest("QED query directory not available")

        mol = lillymol.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O aspirin")
        initial_smiles = mol.smiles()

        calc = lillymol.QED(initialise_from_environment=False)
        self.assertTrue(calc.initialise_from_directory(query_dir))
        value = calc.qed(mol)
        self.assertIsNotNone(value)
        self.assertGreaterEqual(value, 0.0)
        self.assertLessEqual(value, 1.0)
        self.assertEqual(mol.smiles(), initial_smiles)

        calc_from_dir = lillymol.QED(query_dir)
        self.assertAlmostEqual(calc_from_dir.score(mol), value, delta=1.0e-6)

    def test_iwdescr_numpy(self):
        try:
            import gc
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not available")

        calc = lillymol.IWDescr()
        names = calc.feature_names()
        self.assertGreater(len(names), 0)
        self.assertEqual(names, calc.names())
        self.assertEqual(len(names), calc.number_descriptors())

        ethanol = lillymol.MolFromSmiles("CCO ethanol")
        values = calc.process(ethanol)
        self.assertEqual(values.dtype, np.dtype("float32"))
        self.assertEqual(values.shape, (len(names),))
        self.assertEqual(values.ndim, 1)

        mols = [
            lillymol.MolFromSmiles("CCO ethanol"),
            lillymol.MolFromSmiles("c1ccccc1 benzene"),
        ]
        matrix = calc.process_list(mols)
        self.assertEqual(matrix.dtype, np.dtype("float32"))
        self.assertEqual(matrix.shape, (2, len(names)))
        name_to_col = {name: i for i, name in enumerate(names)}
        self.assertEqual(matrix[0, name_to_col["nrings"]], 0.0)
        self.assertEqual(matrix[1, name_to_col["nrings"]], 1.0)

        empty = calc.process_list([])
        self.assertEqual(empty.dtype, np.dtype("float32"))
        self.assertEqual(empty.shape, (0, len(names)))

        del values, matrix, empty
        gc.collect()

    def test_molecular_descriptors_compatibility(self):
        try:
            import gc
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not available")

        calc = lillymol.MolecularDescriptors()
        names = calc.names()
        self.assertGreater(len(names), 0)
        self.assertEqual(names, calc.feature_names())

        mol = lillymol.MolFromSmiles("CCO ethanol")
        array = calc.compute_array(mol)
        self.assertEqual(array.dtype, np.dtype("float32"))
        self.assertEqual(array.shape, (len(names),))

        values = calc.compute(lillymol.MolFromSmiles("CCO ethanol"))
        self.assertEqual(set(values), set(names))
        self.assertEqual(values["nrings"], 0.0)

        matrix = calc.compute_list([
            lillymol.MolFromSmiles("CCO ethanol"),
            lillymol.MolFromSmiles("c1ccccc1 benzene"),
        ])
        self.assertEqual(matrix.dtype, np.dtype("float32"))
        self.assertEqual(matrix.shape, (2, len(names)))
        self.assertEqual(matrix[1, names.index("nrings")], 1.0)

        del array, matrix
        gc.collect()

    def test_jwcats(self):
        charges = _charges_query_dir()
        hbonds = _hbonds_query_dir()
        if charges is None or hbonds is None:
            self.skipTest("JWCats query directories not available")

        calc = lillymol.JWCats(initialise_default_assigners=False)
        self.assertTrue(calc.build_assigners(charges, hbonds))
        self.assertTrue(calc.initialise())
        self.assertTrue(calc.initialised())

        names = calc.feature_names()
        self.assertGreater(len(names), 0)
        mol = lillymol.MolFromSmiles("CCN(CC)C tertiary_amine")
        try:
            import numpy as np
        except ImportError:
            self.skipTest("numpy is not available")

        values = calc.process(mol)
        self.assertEqual(values.dtype, np.dtype("float64"))
        self.assertEqual(values.shape, (len(names),))

        matrix = calc.process_list([
            lillymol.MolFromSmiles("CCN(CC)C tertiary_amine"),
            lillymol.MolFromSmiles("CCO ethanol"),
        ])
        self.assertEqual(matrix.dtype, np.dtype("float64"))
        self.assertEqual(matrix.shape, (2, len(names)))

        calc.set_include_hydrophobic_pairs(False)
        self.assertTrue(calc.initialise())
        self.assertLess(len(calc.feature_names()), len(names))

    def test_rotatable_bonds(self):
        rotbond_calc = lillymol.RotatableBonds()
        rotbond_calc.set_calculation_type(lillymol.EXPENSIVE)

        mol = lillymol.MolFromSmiles("CC ethane")
        self.assertEqual(rotbond_calc.rotatable_bonds(mol), 0)
        mol = lillymol.MolFromSmiles("CCC propane")
        self.assertEqual(rotbond_calc.rotatable_bonds(mol), 0)
        mol = lillymol.MolFromSmiles("CC(F)(F)F trifluoroethane")
        self.assertEqual(rotbond_calc.rotatable_bonds(mol), 0)
        mol = lillymol.MolFromSmiles("CCCC butane")
        self.assertEqual(rotbond_calc.rotatable_bonds(mol), 1)
        self.assertEqual(rotbond_calc.rotatable_bond_atoms(mol), [(1, 2)])
        mol = lillymol.MolFromSmiles("C1CC1C methylcyclopropane")
        self.assertEqual(rotbond_calc.rotatable_bonds(mol), 0)
        self.assertEqual(rotbond_calc.rotatable_bond_atoms(mol), [])

    def test_rotatable_bond_calculation_type(self):
        mol = lillymol.MolFromSmiles("CC(=O)NCC amide")
        rotbond_calc = lillymol.RotatableBonds()
        rotbond_calc.set_calculation_type(lillymol.RotBond.QUICK)
        self.assertEqual(rotbond_calc.rotatable_bond_atoms(mol), [(1, 3), (3, 4)])
        rotbond_calc.set_calculation_type(lillymol.RotBond.EXPENSIVE)
        self.assertEqual(rotbond_calc.rotatable_bond_atoms(mol), [(3, 4)])

    def test_descriptor_helpers(self):
        mol = lillymol.MolFromSmiles("CCO ethanol")
        self.assertIsNotNone(lillymol.alogp(mol))
        self.assertIsNotNone(lillymol.xlogp(mol))
        self.assertIsNotNone(lillymol.tpsa(mol))
        self.assertEqual(lillymol.HbaHbd(mol), (1, 1))
        self.assertEqual(lillymol.NumHAcceptors(mol), mol.lipinski_num_h_acceptors())
        self.assertEqual(lillymol.NumHDonors(mol), mol.lipinski_num_h_donors())
        self.assertEqual(lillymol.RDKitNumHAcceptors(mol), mol.rdkit_num_h_acceptors())
        self.assertEqual(lillymol.RDKitNumHDonors(mol), mol.rdkit_num_h_donors())
        self.assertAlmostEqual(lillymol.fraction_csp3(mol), 1.0)

        benzene = lillymol.MolFromSmiles("c1ccccc1 benzene")
        self.assertAlmostEqual(lillymol.fraction_csp3(benzene), 0.0)
        water = lillymol.MolFromSmiles("O water")
        self.assertAlmostEqual(lillymol.fraction_csp3(water), 0.0)


class TestNanobindTSubstructure(LillyMolNanobindTestCase):

    def test_no_queries(self):
        ts = lillymol.TSubstructure()
        self.assertFalse(ts.substructure_search("C methane"))

    def test_single_query_single_molecule(self):
        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        mol = lillymol.MolFromSmiles("CC")
        self.assertTrue(ts.substructure_search(mol))
        self.assertEqual(ts.num_matches(mol), [2])
        self.assertEqual(ts.number_queries(), 1)
        self.assertEqual(ts.query_names(), ["carbon"])

    def test_multiple_query_single_molecule(self):
        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        mol = lillymol.MolFromSmiles("CC")
        self.assertEqual(ts.num_matches(mol), [2, 0])

    def test_batch_substructure_search_molecules(self):
        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        smiles = ["C methane", "CC ethane", "CCC propane", "C1CC1 cyclopropane", "c1ccccc1 benzene"]
        mols = [lillymol.MolFromSmiles(smi) for smi in smiles]
        self.assertEqual(ts.substructure_search(mols), [True, True, True, True, False])

    def test_batch_substructure_search_smiles(self):
        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        smiles = ["C methane", "N nitrogen", "CN methylamine"]
        self.assertEqual(ts.substructure_search(smiles), [False, True, True])

    def test_batch_num_matches_molecules(self):
        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        smiles = ["C methane", "CC ethane", "N nitrogen", "O oxygen", "CN CN"]
        mols = [lillymol.MolFromSmiles(smi) for smi in smiles]
        self.assertEqual(ts.num_matches(mols), [[1, 0], [2, 0], [0, 1], [0, 0], [1, 1]])

    def test_must_match_all_queries(self):
        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C carbon"))
        self.assertTrue(ts.add_query_from_smarts("N nitrogen"))
        mol = lillymol.MolFromSmiles("C")
        self.assertTrue(ts.substructure_search(mol))
        ts.must_match_all_queries = True
        self.assertFalse(ts.substructure_search(mol))

    def test_aromatic_bonds_kekule_identity(self):
        # The mode is process wide, so put back whatever it was.
        original = lillymol.aromatic_bonds_lose_kekule_identity()
        self.addCleanup(lillymol.set_aromatic_bonds_lose_kekule_identity, original)

        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("a-a biaryl"))
        benzene = lillymol.MolFromSmiles("c1ccccc1 benzene")
        biphenyl = lillymol.MolFromSmiles("c1ccccc1-c1ccccc1 biphenyl")

        # 0: aromatic bonds still match their Kekule single bonds, so a-a is in benzene.
        lillymol.set_aromatic_bonds_lose_kekule_identity(0)
        self.assertEqual(lillymol.aromatic_bonds_lose_kekule_identity(), 0)
        self.assertEqual(ts.substructure_search([benzene, biphenyl]), [True, True])

        # 1 and 2: a-a is only a single bond between two aromatic atoms.
        for mode in (1, 2):
            lillymol.set_aromatic_bonds_lose_kekule_identity(mode)
            self.assertEqual(lillymol.aromatic_bonds_lose_kekule_identity(), mode)
            self.assertEqual(ts.substructure_search([benzene, biphenyl]), [False, True])

        # The non ring single bond form means the same in every mode.
        ts2 = lillymol.TSubstructure()
        self.assertTrue(ts2.add_query_from_smarts("a-!@a biaryl"))
        for mode in (0, 1, 2):
            lillymol.set_aromatic_bonds_lose_kekule_identity(mode)
            self.assertEqual(ts2.substructure_search([benzene, biphenyl]), [False, True])

    def test_aromatic_bonds_kekule_identity_invalid_mode(self):
        original = lillymol.aromatic_bonds_lose_kekule_identity()
        self.addCleanup(lillymol.set_aromatic_bonds_lose_kekule_identity, original)
        for mode in (-1, 3):
            with self.assertRaises(ValueError):
                lillymol.set_aromatic_bonds_lose_kekule_identity(mode)
        self.assertEqual(lillymol.aromatic_bonds_lose_kekule_identity(), original)

    def test_unique_embeddings_only(self):
        ts = lillymol.TSubstructure()
        mol = lillymol.MolFromSmiles("CC(C)(C)c1c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c1C(C)(C)C")
        self.assertTrue(ts.add_query_from_smarts("CC(C)(C)c1c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c1C(C)(C)C"))
        self.assertEqual(ts.num_matches(mol), [559872])
        ts.set_unique_embeddings_only(True)
        self.assertEqual(ts.num_matches(mol), [1])

    def test_max_matches_to_find(self):
        ts = lillymol.TSubstructure()
        mol = lillymol.MolFromSmiles("c1ccccc1")
        self.assertTrue(ts.add_query_from_smarts("c1ccccc1"))
        self.assertEqual(ts.num_matches(mol), [12])
        ts.set_max_matches_to_find(5)
        self.assertEqual(ts.num_matches(mol), [5])

    def test_reduce_to_largest_fragment(self):
        ts = lillymol.TSubstructure()
        mol = lillymol.MolFromSmiles("CC.C")
        self.assertTrue(ts.add_query_from_smarts("C"))
        self.assertEqual(ts.num_matches(mol), [3])
        ts.set_reduce_to_largest_fragment(True)
        self.assertEqual(ts.num_matches(mol), [2])

    def test_make_implicit_hydrogens_explicit(self):
        ts = lillymol.TSubstructure()
        mol = lillymol.MolFromSmiles("CC")
        self.assertTrue(ts.add_query_from_smarts("[#1]-C"))
        self.assertEqual(ts.num_matches(mol), [0])
        ts.set_make_implicit_hydrogens_explicit(True)
        self.assertEqual(ts.num_matches(mol), [6])

    def test_label_matched_atoms(self):
        ts = lillymol.TSubstructure()
        mol = lillymol.MolFromSmiles("Cc1ccccc1")
        self.assertTrue(ts.add_query_from_smarts("C"))
        self.assertTrue(ts.add_query_from_smarts("c"))
        ts.isotope = 4
        self.assertEqual(ts.label_matched_atoms(mol), 2)
        self.assertEqual(mol.unique_smiles(), "[4CH3][4c]1[4cH][4cH][4cH][4cH][4cH]1")

    def test_matched_atoms(self):
        ts = lillymol.TSubstructure()
        self.assertTrue(ts.add_query_from_smarts("C"))
        self.assertTrue(ts.add_query_from_smarts("N"))
        mol = lillymol.MolFromSmiles("CC")
        self.assertEqual(ts.matched_atoms(mol), [[0, 1], []])


if __name__ == "__main__":
    unittest.main()
