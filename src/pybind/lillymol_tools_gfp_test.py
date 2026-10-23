import os
import struct
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from lillymol import Molecule
from lillymol_tools import (
    GFP, GFPContext, GFPList, TruncatedDistanceMatrix,
    TruncatedDistanceMatrixProto, TruncatedDistanceMatrixStorage,
)
from Utilities.GFP_Tools import nearneighbours_pb2


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
    fd, fname = tempfile.mkstemp(suffix='.nn.tfrecord')
    with os.fdopen(fd, 'wb') as output:
        for proto in protos:
            data = proto.SerializeToString()
            length = struct.pack('<Q', len(data))
            output.write(length)
            output.write(struct.pack('<I', _masked_crc32c(length)))
            output.write(data)
            output.write(struct.pack('<I', _masked_crc32c(data)))
    return fname


def _nearneighbours(name, nbrs):
    proto = nearneighbours_pb2.NearNeighbours(name=name)
    for nbr_id, dist in nbrs:
        nbr = proto.nbr.add()
        nbr.id = nbr_id
        nbr.dist = dist
    return proto


def _write_file(contents):
    fd, fname = tempfile.mkstemp(suffix='.gfp')
    with os.fdopen(fd, 'w') as output:
        output.write(contents)
    return fname


def _runfile(path):
    candidates = [
        os.path.join(os.path.dirname(__file__), path),
        os.path.join(os.getcwd(), path),
        os.path.join(os.getcwd(), 'bazel-bin', path),
    ]
    test_srcdir = os.environ.get('TEST_SRCDIR')
    if test_srcdir:
        candidates.extend([
            os.path.join(test_srcdir, '_main', path),
            os.path.join(test_srcdir, path),
        ])
    runfiles_dir = os.environ.get('RUNFILES_DIR')
    if runfiles_dir:
        candidates.extend([
            os.path.join(runfiles_dir, '_main', path),
            os.path.join(runfiles_dir, path),
        ])
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(path)

def _testdata_file(fname):
    candidates = [
        os.path.join(os.path.dirname(__file__), 'testdata', fname),
    ]
    test_srcdir = os.environ.get('TEST_SRCDIR')
    if test_srcdir:
        candidates.extend([
            os.path.join(test_srcdir, '_main', 'pybind', 'testdata', fname),
            os.path.join(test_srcdir, 'pybind', 'testdata', fname),
        ])
    runfiles_dir = os.environ.get('RUNFILES_DIR')
    if runfiles_dir:
        candidates.extend([
            os.path.join(runfiles_dir, '_main', 'pybind', 'testdata', fname),
            os.path.join(runfiles_dir, 'pybind', 'testdata', fname),
        ])
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(fname)

def _mol(smiles):
    mol = Molecule()
    if not mol.build_from_smiles(smiles):
        raise ValueError(smiles)
    return mol


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


class TestTruncatedDistanceMatrix(unittest.TestCase):

    def _run_storage_case(self, storage):
        fname = _write_tfdatarecord([
            _nearneighbours('A', [('B', 0.25), ('C', 0.50)]),
            _nearneighbours('B', [('A', 0.25)]),
            _nearneighbours('C', [('A', 0.50)]),
            _nearneighbours('D', []),
        ])
        try:
            dm = TruncatedDistanceMatrix(fname, storage=storage)
            self.assertEqual(len(dm), 4)
            self.assertEqual(dm.size(), 4)
            self.assertEqual(dm.number_distances(), 2)
            self.assertEqual(dm.index('A'), 0)
            self.assertIsNone(dm.index('missing'))
            self.assertEqual(dm.name(2), 'C')

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
            self.assertEqual(len(dm.distances_or_default([0, 1, 3], [1, 2, 3])), 3)
        finally:
            os.remove(fname)

    def test_row_sparse(self):
        self._run_storage_case(TruncatedDistanceMatrixStorage.ROW_SPARSE)

    def test_row_hash(self):
        self._run_storage_case(TruncatedDistanceMatrixStorage.ROW_HASH)

    def test_one_byte_duplicate_difference_stores_smaller_value(self):
        fname = _write_tfdatarecord([
            _nearneighbours('A', [('B', 10.0 / 255.0)]),
            _nearneighbours('B', [('A', 11.0 / 255.0)]),
        ])
        try:
            dm = TruncatedDistanceMatrix(fname)
            self.assertEqual(dm.number_distances(), 1)
            self.assertEqual(dm.duplicate_distances_differing_by_one(), 1)
            self.assertAlmostEqual(dm.distance(0, 1), 10.0 / 255.0, places=6)
        finally:
            os.remove(fname)

    def test_tbb_tool_writes_indexed_tfdatarecord(self):
        exe = _runfile('Utilities/GFP_Tools/gfp_nearneighbours_single_file_tbb')
        input_fname = _testdata_file('rand10.standard.gfp')
        fd, output_fname = tempfile.mkstemp(suffix='.indexed.nn.tfdata')
        os.close(fd)
        os.remove(output_fname)
        try:
            subprocess.run([
                exe, '-T', '0.6', '-S', 'indexed', '-S', output_fname, input_fname,
            ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            dm = TruncatedDistanceMatrix(
                output_fname,
                storage=TruncatedDistanceMatrixStorage.ROW_SPARSE,
                proto_type=TruncatedDistanceMatrixProto.NEARNEIGHBOURS_INDICES,
            )
            self.assertEqual(dm.size(), 12)
            self.assertGreater(dm.number_distances(), 0)
            self.assertEqual(dm.index(dm.name(0)), 0)
            self.assertAlmostEqual(dm.distance_or_default(0, 0), 0.0, places=6)
        finally:
            if os.path.exists(output_fname):
                os.remove(output_fname)

    def test_rejects_conflicting_duplicate_distance(self):
        fname = _write_tfdatarecord([
            _nearneighbours('A', [('B', 0.25)]),
            _nearneighbours('B', [('A', 0.75)]),
        ])
        try:
            with self.assertRaises(RuntimeError):
                TruncatedDistanceMatrix(fname)
        finally:
            os.remove(fname)


class TestGFPList(unittest.TestCase):

    def test_read_and_distance(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            self.assertEqual(len(gfp), 2)
            self.assertEqual(gfp.size(), 2)
            self.assertEqual(gfp.tags(), ['FCTS<'])
            self.assertEqual(gfp.smiles(0), 'C')
            self.assertEqual(gfp.id(1), 'Ethane')
            self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
            self.assertAlmostEqual(gfp.distance(1, 0), 0.5, places=6)
        finally:
            os.remove(fname)

    def test_nearest_neighbours(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            hits = gfp.nearest_neighbours(0, 1)
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0].index, 1)
            self.assertAlmostEqual(hits[0].distance, 0.5, places=6)

            close = gfp.nearest_neighbours_within_distance(0, 0.5)
            self.assertEqual(len(close), 1)
            self.assertEqual(close[0].index, 1)
            self.assertAlmostEqual(close[0].distance, 0.5, places=6)

            none = gfp.nearest_neighbours_within_distance(0, 0.49)
            self.assertEqual(none, [])
        finally:
            os.remove(fname)

    def test_standard_gfp_golden_distances(self):
        gfp = GFPList.from_file(_testdata_file('rand10.standard.gfp'))
        self.assertGreater(len(gfp), 3)
        self.assertEqual(gfp.id(0), 'CHEMBL3460651')
        self.assertEqual(gfp.id(1), 'CHEMBL3460651.a')
        self.assertEqual(gfp.id(3), 'CHEMBL1417367')

        self.assertAlmostEqual(gfp.distance(0, 1), 0.0421, delta=0.0001)
        self.assertAlmostEqual(gfp.distance(1, 0), 0.0421, delta=0.0001)
        self.assertAlmostEqual(gfp.distance(3, 0), 0.499, delta=0.001)
        self.assertAlmostEqual(gfp.distance(0, 3), 0.499, delta=0.001)

    def test_complex_gfp_file(self):
        gfp = GFPList.from_file(_testdata_file('rand10.complex.gfp'))
        self.assertEqual(len(gfp), 10)
        tags = gfp.tags()
        self.assertGreaterEqual(len(tags), 9)
        self.assertIn('FPIW<', tags)
        self.assertIn('FPMK<', tags)
        self.assertIn('MPR<', tags)

        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
        d01 = gfp.distance(0, 1)
        self.assertGreaterEqual(d01, 0.0)
        self.assertLessEqual(d01, 1.0)
        self.assertAlmostEqual(d01, gfp.distance(1, 0), places=6)

        hits = gfp.nearest_neighbours(0, 3)
        self.assertEqual(len(hits), 3)
        self.assertNotIn(0, [hit.index for hit in hits])
        self.assertEqual([hit.distance for hit in hits], sorted(hit.distance for hit in hits))

        within = gfp.nearest_neighbours_within_distance(0, 1.0)
        self.assertEqual(len(within), 9)
        self.assertNotIn(0, [hit.index for hit in within])

    def test_use_only_and_use_all(self):
        fname = _write_file(CARBON_GFP)
        try:
            gfp = GFPList.from_file(fname)
            gfp.use_only(['FCTS<'])
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
            gfp.use_all()
            self.assertAlmostEqual(gfp.distance(0, 1), 0.5, places=6)
        finally:
            os.remove(fname)

    def test_generator_specs_match_standard_context(self):
        standard = GFPContext.standard()
        from_specs = GFPContext.from_specs([GFP.iw(), GFP.maccs(), GFP.mpr()])

        self.assertEqual(from_specs.tags(), standard.tags())

        mol = _mol('CCO ethanol')
        fp_standard = standard.fingerprint(mol)
        fp_from_specs = from_specs.fingerprint(mol)
        self.assertAlmostEqual(standard.distance(fp_standard, fp_from_specs),
                               0.0, places=6)

    def test_generator_specs_are_order_independent(self):
        ctx1 = GFPContext.from_specs([GFP.iw(), GFP.maccs(), GFP.mpr()])
        ctx2 = GFPContext.from_specs([GFP.mpr(), GFP.maccs(), GFP.iw()])

        self.assertEqual(ctx1.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])
        self.assertEqual(ctx1.tags(), ctx2.tags())

        mol = _mol('CCC propane')
        fp1 = ctx1.fingerprint(mol)
        fp2 = ctx2.fingerprint(mol)
        self.assertAlmostEqual(ctx1.distance(fp1, fp2), 0.0, places=6)

    def test_generator_specs_reject_duplicates(self):
        with self.assertRaises(RuntimeError):
            GFPContext.from_specs([GFP.iw(), GFP.iw()])

    def test_maccs_level2_false(self):
        spec = GFP.maccs(level2=False)
        self.assertEqual(spec.components(), ['FPMK<'])
        self.assertEqual(repr(spec), 'GFP.maccs(level2=False)')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['FPMK<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)





    def test_xlogp_generator_spec(self):
        spec = GFP.xlogp(replicates=4)
        self.assertEqual(spec.components(), ['NCXLOGP4<'])
        self.assertEqual(repr(spec), 'GFP.xlogp(replicates=4)')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCXLOGP4<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_xlogp_rejects_invalid_replicates(self):
        with self.assertRaises(ValueError):
            GFP.xlogp(replicates=0)

    def test_tpsa_generator_spec(self):
        spec = GFP.tpsa(replicates=4)
        self.assertEqual(spec.components(), ['NCTPSA4<'])
        self.assertEqual(repr(spec), 'GFP.tpsa(replicates=4)')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCTPSA4<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_tpsa_rejects_invalid_replicates(self):
        with self.assertRaises(ValueError):
            GFP.tpsa(replicates=0)

    def test_formula_generator_spec(self):
        spec = GFP.formula()
        self.assertEqual(spec.components(), ['FCFML<'])
        self.assertEqual(repr(spec), 'GFP.formula()')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['FCFML<'])
        ethane = _mol('CC ethane')
        propane = _mol('CCC propane')
        fp1 = ctx.fingerprint(ethane)
        fp2 = ctx.fingerprint(propane)
        self.assertAlmostEqual(ctx.distance(fp1, fp1), 0.0, places=6)
        self.assertGreater(ctx.distance(fp1, fp2), 0.0)
        self.assertLess(ctx.distance(fp1, fp2), 1.0)

    def test_cats_generator_spec(self):
        spec = GFP.cats(max_path_length=4)
        self.assertEqual(spec.components(), ['NCCATS4<'])
        self.assertEqual(
            repr(spec),
            'GFP.cats(max_path_length=4, include_hydrophobic_pairs=True)')

        if not os.environ.get('LILLYMOL_HOME'):
            self.skipTest('LILLYMOL_HOME not set')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCCATS4<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_catsp_generator_spec(self):
        spec = GFP.cats(max_path_length=4, include_hydrophobic_pairs=False)
        self.assertEqual(spec.components(), ['NCCATSP4<'])
        self.assertEqual(
            repr(spec),
            'GFP.cats(max_path_length=4, include_hydrophobic_pairs=False)')

        if not os.environ.get('LILLYMOL_HOME'):
            self.skipTest('LILLYMOL_HOME not set')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCCATSP4<'])

    def test_cats_rejects_invalid_max_path_length(self):
        with self.assertRaises(ValueError):
            GFP.cats(max_path_length=0)

    def test_atom_pair_generator_spec(self):
        spec = GFP.atom_pair(min_separation=1, max_separation=3, atom_type='UST:Y')
        self.assertEqual(spec.components(), ['NCAP1M3USTY<'])
        self.assertEqual(
            repr(spec),
            "GFP.atom_pair(min_separation=1, max_separation=3, "
            "atom_type='UST:Y', include_out_of_range=False)")

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCAP1M3USTY<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_atom_pair_out_of_range_tag(self):
        spec = GFP.atom_pair(
            min_separation=0, max_separation=2, include_out_of_range=True)
        self.assertEqual(spec.components(), ['NCAPT0M2USTY<'])

    def test_atom_pair_rejects_invalid_specs(self):
        with self.assertRaises(ValueError):
            GFP.atom_pair(min_separation=-1)
        with self.assertRaises(ValueError):
            GFP.atom_pair(min_separation=4, max_separation=3)
        with self.assertRaises(ValueError):
            GFP.atom_pair(atom_type='')
        with self.assertRaises(RuntimeError):
            GFPContext.from_specs([GFP.atom_pair(atom_type='UST:Y-')])

    def test_ring_substitution_generator_spec(self):
        spec = GFP.ring_substitution()
        self.assertEqual(spec.components(), ['NCRS<'])
        self.assertEqual(repr(spec), 'GFP.ring_substitution()')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCRS<'])
        mol = _mol('c1ccccc1C toluene')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_spinach_generator_spec(self):
        spec = GFP.spinach()
        self.assertEqual(spec.components(), ['FPSPIN<'])
        self.assertEqual(repr(spec), 'GFP.spinach(label_join_points=False)')

        labelled = GFP.spinach(label_join_points=True)
        self.assertEqual(labelled.components(), ['FPSPINI<'])
        self.assertEqual(repr(labelled), 'GFP.spinach(label_join_points=True)')

        ctx = GFPContext.from_specs([spec, labelled])
        self.assertEqual(ctx.tags(), ['FPSPIN<', 'FPSPINI<'])
        mol = _mol('c1ccccc1C toluene')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_scaffold_generator_spec(self):
        spec = GFP.scaffold()
        self.assertEqual(spec.components(), ['FPSCAF<'])
        self.assertEqual(repr(spec), 'GFP.scaffold(label_join_points=False)')

        labelled = GFP.scaffold(label_join_points=True)
        self.assertEqual(labelled.components(), ['FPSCAFI<'])
        self.assertEqual(repr(labelled), 'GFP.scaffold(label_join_points=True)')

        ctx = GFPContext.from_specs([spec, labelled])
        self.assertEqual(ctx.tags(), ['FPSCAF<', 'FPSCAFI<'])
        mol = _mol('c1ccccc1C toluene')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_scaffold_and_spinach_together(self):
        ctx = GFPContext.from_specs([
            GFP.spinach(),
            GFP.spinach(label_join_points=True),
            GFP.scaffold(),
            GFP.scaffold(label_join_points=True),
        ])
        self.assertEqual(ctx.tags(), ['FPSCAF<', 'FPSCAFI<', 'FPSPIN<', 'FPSPINI<'])
        mol = _mol('c1ccccc1C toluene')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_substructure_generator_spec(self):
        spec = GFP.substructure('c1ccccc1', radius=1)
        components = spec.components()
        self.assertEqual(len(components), 1)
        self.assertTrue(components[0].startswith('FPSUB1USTARY'))
        self.assertTrue(components[0].endswith('<'))
        self.assertEqual(
            repr(spec),
            "GFP.substructure(smarts='c1ccccc1', radius=1, "
            "atom_type='UST:ARY', no_match='empty')")

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), components)
        benzene = _mol('c1ccccc1 benzene')
        fp = ctx.fingerprint(benzene)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_substructure_radius_changes_tag(self):
        radius0 = GFP.substructure('c1ccccc1', radius=0)
        radius2 = GFP.substructure('c1ccccc1', radius=2)
        self.assertNotEqual(radius0.components(), radius2.components())
        ctx = GFPContext.from_specs([radius0, radius2])
        self.assertEqual(len(ctx.tags()), 2)

    def test_substructure_no_match_empty(self):
        ctx = GFPContext.from_specs([GFP.substructure('N', no_match='empty')])
        methane = _mol('C methane')
        fp = ctx.fingerprint(methane)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_substructure_no_match_error(self):
        ctx = GFPContext.from_specs([GFP.substructure('N', no_match='error')])
        methane = _mol('C methane')
        with self.assertRaises(RuntimeError):
            ctx.fingerprint(methane)

    def test_substructure_rejects_invalid_specs(self):
        with self.assertRaises(ValueError):
            GFP.substructure('')
        with self.assertRaises(ValueError):
            GFP.substructure('C', radius=-1)
        with self.assertRaises(ValueError):
            GFP.substructure('C', atom_type='')
        with self.assertRaises(ValueError):
            GFP.substructure('C', no_match='skip')
        with self.assertRaises(RuntimeError):
            GFPContext.from_specs([GFP.substructure('[')])
        with self.assertRaises(RuntimeError):
            GFPContext.from_specs([GFP.substructure('C', atom_type='UST:A-Y')])

    def test_ec_generator_spec(self):
        spec = GFP.ec(radius=3, atom_type='UST:AY')
        self.assertEqual(spec.components(), ['NCEC3USTAY<'])
        self.assertEqual(repr(spec), "GFP.ec(radius=3, atom_type='UST:AY')")

        ctx = GFPContext.from_specs([GFP.ec(), spec])
        self.assertEqual(ctx.tags(), ['NCEC3USTAY<', 'NCEC3USTZ<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_ec_rejects_invalid_specs(self):
        with self.assertRaises(ValueError):
            GFP.ec(radius=-1)
        with self.assertRaises(ValueError):
            GFP.ec(atom_type='')
        with self.assertRaises(RuntimeError):
            GFPContext.from_specs([GFP.ec(radius=3, atom_type='UST:A-Y')])

    def test_alogp_generator_spec(self):
        spec = GFP.alogp(replicates=4)
        self.assertEqual(spec.components(), ['NCALOGP4<'])
        self.assertEqual(repr(spec), 'GFP.alogp(replicates=4)')

        ctx = GFPContext.from_specs([spec])
        self.assertEqual(ctx.tags(), ['NCALOGP4<'])
        mol = _mol('CCO ethanol')
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

    def test_alogp_rejects_invalid_replicates(self):
        with self.assertRaises(ValueError):
            GFP.alogp(replicates=0)

    def test_standard_context_generates_fingerprints(self):
        ctx = GFPContext.standard()
        self.assertEqual(ctx.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])
        self.assertTrue(ctx.can_generate_fingerprints())

        ethane = Molecule()
        self.assertTrue(ethane.build_from_smiles('CC ethane'))
        propane = Molecule()
        self.assertTrue(propane.build_from_smiles('CCC propane'))

        fp1 = ctx.fingerprint(ethane)
        fp2 = ctx.fingerprint(propane)
        self.assertAlmostEqual(ctx.distance(fp1, fp1), 0.0, places=6)
        d12 = ctx.distance(fp1, fp2)
        self.assertGreaterEqual(d12, 0.0)
        self.assertLessEqual(d12, 1.0)
        self.assertAlmostEqual(ctx.distance(fp2, fp1), d12, places=6)

    def test_standard_list_add_and_query_fingerprint(self):
        gfp = GFPList.standard()
        self.assertEqual(gfp.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])

        for smiles in ['CC ethane', 'CCC propane', 'CCCC butane']:
            mol = Molecule()
            self.assertTrue(mol.build_from_smiles(smiles))
            gfp.add(mol)

        self.assertEqual(len(gfp), 3)
        self.assertEqual(gfp.id(0), 'ethane')
        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)

        query = Molecule()
        self.assertTrue(query.build_from_smiles('CCC query'))
        fp = GFPContext.standard().fingerprint(query)
        self.assertAlmostEqual(gfp.distance(fp, 1), 0.0, places=6)

        hits = gfp.nearest_neighbours(fp, 2)
        self.assertEqual(hits[0].index, 1)
        self.assertAlmostEqual(hits[0].distance, 0.0, places=6)


    def test_programmatic_standard_list_matches_stored_file_distances(self):
        stored = GFPList.from_file(_testdata_file('rand10.standard.gfp'))
        molecules = [_mol(stored.smiles(i)) for i in range(len(stored))]
        generated = GFPList.standard_from_molecules(molecules)

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
                                       msg=f'{i},{j}')

    def test_standard_from_molecules_can_store_metadata(self):
        molecules = [_mol('CC ethane'), _mol('CCC propane'), _mol('CCCC butane')]
        gfp = GFPList.standard_from_molecules(molecules, store_metadata=True)

        self.assertEqual(len(gfp), 3)
        self.assertTrue(gfp.metadata_stored())
        self.assertEqual(gfp.id(0), 'ethane')
        self.assertEqual(gfp.smiles(1), 'CCC')
        self.assertAlmostEqual(gfp.distance(1, 1), 0.0, places=6)

    def test_add_molecules_defaults_to_no_metadata(self):
        molecules = [_mol('CC ethane'), _mol('CCC propane')]
        gfp = GFPList.standard()
        gfp.add_molecules(molecules)

        self.assertEqual(len(gfp), 2)
        self.assertFalse(gfp.metadata_stored())
        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)
        with self.assertRaises(RuntimeError):
            gfp.id(0)
        with self.assertRaises(RuntimeError):
            gfp.add(_mol('CCCC butane'))

    def test_gfp_failures_raise_exceptions(self):
        mol = Molecule()
        self.assertTrue(mol.build_from_smiles('CC ethane'))

        with self.assertRaises(RuntimeError):
            GFPContext().fingerprint(mol)

        with self.assertRaises(RuntimeError):
            GFPList().add(mol)

        with self.assertRaises(RuntimeError):
            GFPList.from_file('/no/such/gfp/file.gfp')

        gfp = GFPList.standard()
        gfp.add(mol)
        with self.assertRaises(IndexError):
            gfp.distance(0, 1)
        with self.assertRaises(IndexError):
            gfp.smiles(1)
        with self.assertRaises(IndexError):
            gfp.nearest_neighbours(1, 1)
        with self.assertRaises(ValueError):
            gfp.nearest_neighbours_within_distance(0, -1.0)
        with self.assertRaises(RuntimeError):
            gfp.use_only(['NO_SUCH_TAG<'])

        carbon_fname = _write_file(CARBON_GFP)
        try:
            carbon = GFPList.from_file(carbon_fname)
            fp = GFPContext.standard().fingerprint(mol)
            with self.assertRaises(ValueError):
                carbon.distance(fp, 0)
            with self.assertRaises(ValueError):
                carbon.nearest_neighbours(fp, 1)
        finally:
            os.remove(carbon_fname)


    def test_generated_standard_fingerprints_match_stored_file(self):
        stored = GFPList.from_file(_testdata_file('rand10.standard.gfp'))
        self.assertEqual(stored.tags(), ['FPIW<', 'FPMK<', 'FPMK2<', 'MPR<'])

        ctx = GFPContext.standard()
        for i in range(len(stored)):
            mol = Molecule()
            self.assertTrue(mol.build_from_smiles(stored.smiles(i)))
            fp = ctx.fingerprint(mol)
            self.assertAlmostEqual(stored.distance(fp, i), 0.0, places=6,
                                   msg=stored.id(i))

    def test_standard_generation_can_skip_preprocessing(self):
        ctx = GFPContext.standard(preprocess=False)
        self.assertTrue(ctx.can_generate_fingerprints())
        mol = Molecule()
        self.assertTrue(mol.build_from_smiles('CCO ethanol'))
        fp = ctx.fingerprint(mol)
        self.assertAlmostEqual(ctx.distance(fp, fp), 0.0, places=6)

        gfp = GFPList.standard(preprocess=False)
        gfp.add(mol)
        self.assertEqual(len(gfp), 1)
        self.assertAlmostEqual(gfp.distance(0, 0), 0.0, places=6)


if __name__ == '__main__':
    unittest.main()
