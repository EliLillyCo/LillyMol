# tests for LillyMol io functionality

# Right now this is sub-optimal. When we get py_test working
# we should be able to specify data files within BUILD

from collections import defaultdict
import os
import tempfile

from absl import app
from absl import logging
from absl.testing import absltest

from google.protobuf import text_format

from lillymol import *
from lillymol_io import *

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

SDF = """methane
iwcorina  09042307323D 1   1.00000     0.00000     0
CORINA-API 3.49 0006  12.02.2015
  1  0  0  0  0  0            999 V2000
   -0.0127    1.0858    0.0080 C   0  0  0
M  END
$$$$
ethane
iwcorina  09042307323D 1   1.00000     0.00000     0
CORINA-API 3.49 0006  12.02.2015
  2  1  0  0  0  0            999 V2000
   -0.0187    1.5258    0.0104 C   0  0  0
    0.0021   -0.0041    0.0020 C   0  0  0
  1  2  1  0
M  END
$$$$
"""

ENAMINE="""
  -ISIS-  -- StrEd -- 

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.7500    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  END
>  <idnumber> (Z33546370)
Z33546370

>  <LogS> (Z33546370)
0.5

>  <LogP> (Z33546370)
-1.114

>  <PSA> (Z33546370)
43.09

>  <link> (Z33546370)
https://www.enaminestore.com/catalog/Z33546370

$$$$
"""

SDF_MULTILINE="""methane
iwcorina  09042307323D 1   1.00000     0.00000     0
CORINA-API 3.49 0006  12.02.2015
  1  0  0  0  0  0            999 V2000
   -0.0127    1.0858    0.0080 C   0  0  0
M  END
>  <COMMENT>
first line
second line

>  <SINGLE>
value

$$$$
"""


class TestLillyMolSubstructure(absltest.TestCase):
  def test_open_file_ok_suffix(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    reader = Reader()
    self.assertTrue(reader.open(fname))

  def test_open_file_bad_suffix(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.foo")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    reader = Reader()
    self.assertFalse(reader.open(fname))

  def test_open_file_bad_suffix_provide_file_type(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.foo")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    reader = Reader()
    self.assertTrue(reader.open(fname, FileType.SMI))

  def test_read_by_iter(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    reader = Reader()
    self.assertTrue(reader.open(fname, FileType.SMI))
    molecules_read = 0
    for mol in reader:
      molecules_read += 1

    self.assertEqual(reader.molecules_read(), molecules_read)

  def test_read_smiles_via_loop(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    reader = Reader()
    self.assertTrue(reader.open(fname))
    while True:
      m = reader.next()
      if m is None:
        break;
    self.assertEqual(reader.molecules_read(), 10)

  def test_read_sdf_via_loop(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(SDF)

    reader = Reader()
    self.assertTrue(reader.open(fname))
    molecules_read = 0
    for mol in reader:
      molecules_read += 1
    self.assertEqual(reader.molecules_read(), molecules_read)

  def test_context_aliases(self):
    self.assertIs(ReaderContext, MolReaderContext)
    self.assertIs(ContextWriter, MolWriterContext)

  def test_read_smiles_via_context_reader(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    molecules_read = 0
    with MolReaderContext(fname, FileType.SMI) as reader:
      for mol in reader:
        molecules_read += 1
    self.assertEqual(molecules_read, 10)

  def test_connection_table_errors(self):
    f = SMILES.split('\n')
    f[4] = "foo bar"
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write('\n'.join(f) + '\n')

    set_display_smiles_interpretation_error_messages(False)
    molecules_read = 0
    with ReaderContext(fname, FileType.SMI) as reader:
      for mol in reader:
        molecules_read += 1
    self.assertEqual(molecules_read, 4)

  def test_ignore_connection_table_errors(self):
    f = SMILES.split('\n')
    f[4] = "foo bar"
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write('\n'.join(f) + '\n')

    set_display_smiles_interpretation_error_messages(0)
    molecules_read = 0
    with ReaderContext(fname, FileType.SMI) as reader:
      reader.set_ignore_connection_table_errors(1)
      for mol in reader:
        molecules_read += 1
    self.assertEqual(molecules_read, 9)
    self.assertEqual(reader.connection_table_errors_encountered(), 1)

  def test_molecules_remaining_sdf(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(SDF)

    reader = Reader()
    self.assertTrue(reader.open(fname))
    self.assertEqual(reader.molecules_remaining(), 2)

  def test_molecules_remaining_smi(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    reader = Reader()
    with ReaderContext(fname, FileType.SMI) as reader:
      for i in range(4):
        mol = reader.next()

      self.assertEqual(reader.molecules_remaining(), 6)

  # Writer objects can simultaneously write multiple types.
  def test_write_smiles(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    stem = os.path.join(tmpdir, "anything")
    writer = Writer()
    writer.add_output_type(FileType.SMI)
    writer.add_output_type(FileType.SDF)
    self.assertTrue(writer.new_stem(stem))
    for i in range(4):
      m = LillyMolFromSmiles("C" * (i + 1))
      writer.write(m)
    writer.close()

    suffix_type = {"smi":FileType.SMI, "sdf":FileType.SDF}

    # Keep track of number of times each structure is read.
    seen = defaultdict(int)

    for suffix,ftype in suffix_type.items():
      fname = f"{stem}.{suffix}"
      # logging.info("Opening %s", fname)
      with ReaderContext(fname, ftype) as reader:
        for mol in reader:
          seen[mol.unique_smiles()] += 1

    # Reading smiles and .sdf should yield same results.
    for smi,count in seen.items():
      self.assertEqual(count, 2)

  def test_context_writer_smi(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    stem = os.path.join(tmpdir, "somefile")

    mols = [LillyMolFromSmiles("CN" * (i + 1)) for i in range(10)]

    with MolWriterContext(stem, FileType.SMI) as writer:
      for m in mols:
        writer.write(m)

    fname = os.path.join(tmpdir, "somefile.smi")
    molecules_read = 0
    with ReaderContext(fname) as reader:
      for mol in reader:
        molecules_read += 1

    self.assertEqual(molecules_read, 10)

  def test_sdfid_with_prepend(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    self.assertTrue(set_sdf_identifier("idnumber"))
    set_mdlquiet(True)
    set_ignore_bad_m(True)
    set_prepend_sdfid(True)
    with ReaderContext(fname) as reader:
      for mol in reader:
        self.assertEqual(mol.name(), "idnumber:Z33546370")

    # Reset back to default value.
    self.assertTrue(set_sdf_identifier(""))

  def test_sdfid_no_prepend(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    self.assertTrue(set_sdf_identifier("idnumber"))
    set_mdlquiet(True)
    set_ignore_bad_m(True)
    set_prepend_sdfid(False)
    with ReaderContext(fname) as reader:
      for mol in reader:
        self.assertEqual(mol.name(), "Z33546370")

    # Reset back to default value.
    self.assertTrue(set_sdf_identifier(""))
    set_prepend_sdfid(True)

  def test_sdfid_allsdfid(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    set_allsdfid(True)
    with ReaderContext(fname) as reader:
      for mol in reader:
        self.assertEqual(mol.name(), "idnumber:Z33546370 LogS:0.5 LogP:-1.114 PSA:43.09 link:https://www.enaminestore.com/catalog/Z33546370")

    # Reset back to default value.
    self.assertTrue(set_sdf_identifier(""))
    set_allsdfid(False)

  def test_reader_context_sdf_identifier_keyword(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname, sdf_identifier="idnumber") as reader:
      mol = next(reader)
      self.assertEqual(mol.name(), "idnumber:Z33546370")

    with ReaderContext(fname) as reader:
      mol = next(reader)
      self.assertEqual(mol.name(), "")

  def test_reader_context_sdf_identifier_keyword_no_prepend(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname, sdf_identifier="idnumber", prepend_sdfid=False) as reader:
      mol = next(reader)
      self.assertEqual(mol.name(), "Z33546370")

  def test_reader_context_sdf_options_restore_previous_global(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    self.assertTrue(set_sdf_identifier("LogS"))
    set_prepend_sdfid(False)
    with ReaderContext(fname, sdf_identifier="idnumber") as reader:
      mol = next(reader)
      self.assertEqual(mol.name(), "idnumber:Z33546370")

    with ReaderContext(fname) as reader:
      mol = next(reader)
      self.assertEqual(mol.name(), "0.5")

    self.assertTrue(set_sdf_identifier(""))
    set_prepend_sdfid(True)

  def test_reader_context_independent_sdf_options(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname, sdf_identifier="idnumber", prepend_sdfid=False) as reader1:
      with ReaderContext(fname, sdf_identifier="LogS", prepend_sdfid=False) as reader2:
        self.assertEqual(next(reader1).name(), "Z33546370")
        self.assertEqual(next(reader2).name(), "0.5")

  def test_reader_context_keep_sdf_tags_keyword(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname, keep_sdf_tags=True) as reader:
      mol = next(reader)

    self.assertEqual(mol.sdf_tags()["idnumber"], "Z33546370")
    self.assertEqual(mol.sdf_tags()["LogP"], "-1.114")
    self.assertEqual(mol.sdf_tags()["PSA"], "43.09")

  def test_reader_context_sdf_tags_default_empty(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname) as reader:
      mol = next(reader)

    self.assertEqual(mol.sdf_tags(), {})

  def test_reader_context_sdf_tags_multiline_values(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(SDF_MULTILINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname, keep_sdf_tags=True) as reader:
      mol = next(reader)

    self.assertEqual(mol.sdf_tags()["COMMENT"], "first line\nsecond line")
    self.assertEqual(mol.sdf_tags()["SINGLE"], "value")

  def test_reader_context_reader_options_sdf_identifier(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    options = ReaderOptions(sdf_identifier="LogP", prepend_sdfid=False)
    with ReaderContext(fname, options=options) as reader:
      self.assertEqual(next(reader).name(), "-1.114")

  def test_reader_context_reader_options_with_file_type(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    options = ReaderOptions(sdf_identifier="PSA", prepend_sdfid=False)
    with ReaderContext(fname, FileType.SDF, options=options) as reader:
      self.assertEqual(next(reader).name(), "43.09")

  def test_reader_options_reused(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    options = ReaderOptions(sdf_identifier="idnumber", prepend_sdfid=False)
    with ReaderContext(fname, options=options) as reader1:
      self.assertEqual(next(reader1).name(), "Z33546370")
    with ReaderContext(fname, options=options) as reader2:
      self.assertEqual(next(reader2).name(), "Z33546370")

  def test_reader_options_mutable_fields(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    options = ReaderOptions()
    options.sdf_identifier = "LogS"
    options.prepend_sdfid = False
    with ReaderContext(fname, options=options) as reader:
      self.assertEqual(next(reader).name(), "0.5")

  def test_reader_options_invalid_sdf_identifier_raises(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    options = ReaderOptions(sdf_identifier="[")
    with self.assertRaises(ValueError):
      ReaderContext(fname, options=options)

    with ReaderContext(fname) as reader:
      with self.assertRaises(ValueError):
        reader.apply_options(options)

  def test_reader_context_sdf_tags_to_json_keywords(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname, sdf_tags_to_json=True, all_sdf_tags=True) as reader:
      mol = next(reader)
      self.assertEqual(mol.name(), '{ "idnumber": "Z33546370", "LogS": "0.5", "LogP": "-1.114", "PSA": "43.09", "link": "https://www.enaminestore.com/catalog/Z33546370" }')

  def test_reader_context_set_sdf_options(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    with ReaderContext(fname) as reader:
      reader.set_sdf_options(sdf_identifier="LogS", prepend_sdfid=False)
      mol = next(reader)
      self.assertEqual(mol.name(), "0.5")

  def test_sdfid_sdf_tags_to_json_nothing(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    set_sdf_tags_to_json(True)
    self.assertTrue(set_sdf_identifier(""))
    with ReaderContext(fname) as reader:
      for mol in reader:
        self.assertEqual(mol.name(), "{ }")

    # Reset back to default value.
    self.assertTrue(set_sdf_identifier(""))
    set_sdf_tags_to_json(False)

  def test_sdfid_sdf_tags_to_json_sdf_identifier(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    set_sdf_tags_to_json(True)
    self.assertTrue(set_sdf_identifier("idnumber"))
    with ReaderContext(fname) as reader:
      for mol in reader:
        self.assertEqual(mol.name(), "{ \"idnumber\": \"Z33546370\" }")

    # Reset back to default value.
    self.assertTrue(set_sdf_identifier(""))
    set_sdf_tags_to_json(False)

  def test_sdfid_sdf_tags_to_json_all_sdf_tags(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    set_sdf_tags_to_json(True)
    self.assertTrue(set_sdf_identifier(""))
    set_allsdfid(True)
    with ReaderContext(fname) as reader:
      for mol in reader:
        self.assertEqual(mol.name(), '{ "idnumber": "Z33546370", "LogS": "0.5", "LogP": "-1.114", "PSA": "43.09", "link": "https://www.enaminestore.com/catalog/Z33546370" }')

    # Reset back to default value.
    self.assertTrue(set_sdf_identifier(""))
    set_sdf_tags_to_json(False)
    set_allsdfid(False)

  def test_firstsdfid(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    set_firstsdftag(True)
    self.assertTrue(set_sdf_identifier(""))
    # This has no effect here. Is that a feature or a bug?
    set_prepend_sdfid(True)
    with ReaderContext(fname) as reader:
      for mol in reader:
        self.assertEqual(mol.name(), 'Z33546370')

    # Reset back to default value.
    self.assertTrue(set_sdf_identifier(""))
    set_allsdfid(False)
    set_firstsdftag(False)


  def test_molecule_preprocessing_process(self):
    prep = MoleculePreprocessing(largest_fragment=True, remove_chirality=True,
                                 remove_isotopes=True)
    self.assertTrue(prep.active())

    mol = LillyMolFromSmiles("[1CH3].[2CH]([C@H](N)F)C")
    self.assertEqual(mol.number_fragments(), 2)
    self.assertGreater(mol.number_chiral_centres(), 0)

    self.assertGreater(prep.process(mol), 0)
    self.assertEqual(mol.number_fragments(), 1)
    self.assertEqual(mol.number_chiral_centres(), 0)
    for atom in range(mol.natoms()):
      self.assertEqual(mol.isotope(atom), 0)

  def test_molecule_preprocessing_process_copy(self):
    prep = MoleculePreprocessing(largest_fragment=True)
    mol = LillyMolFromSmiles("C.CC")

    copy = prep.process_copy(mol)
    self.assertEqual(mol.number_fragments(), 2)
    self.assertEqual(copy.number_fragments(), 1)
    self.assertEqual(copy.natoms(), 2)

  def test_reader_context_preprocessing_keywords(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write("C.CC frag\n")

    with ReaderContext(fname, FileType.SMI, largest_fragment=True) as reader:
      mol = next(reader)
      self.assertEqual(mol.number_fragments(), 1)
      self.assertEqual(mol.natoms(), 2)
      self.assertIsNone(reader.next())
      self.assertEqual(reader.molecules_read(), 1)

  def test_reader_options_keep_sdf_tags(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.sdf")
    with open(fname, "w") as writer:
      writer.write(ENAMINE)

    set_mdlquiet(True)
    set_ignore_bad_m(True)
    options = ReaderOptions(keep_sdf_tags=True)
    with ReaderContext(fname, options=options) as reader:
      mol = next(reader)

    self.assertEqual(mol.sdf_tags()["idnumber"], "Z33546370")

  def test_reader_options_preprocessing(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write("C.CC frag\n")

    options = ReaderOptions(largest_fragment=True)
    with ReaderContext(fname, options=options) as reader:
      mol = next(reader)
      self.assertEqual(mol.number_fragments(), 1)
      self.assertEqual(mol.natoms(), 2)

  def test_reader_context_set_preprocessing(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write("[1CH4] iso\n")

    with ReaderContext(fname) as reader:
      self.assertFalse(reader.preprocessing_active())
      reader.set_preprocessing(remove_isotopes=True)
      self.assertTrue(reader.preprocessing_active())
      mol = next(reader)
      self.assertEqual(mol.isotope(0), 0)

if __name__ == '__main__':
  #app.run(absltest.main)
  absltest.main()
