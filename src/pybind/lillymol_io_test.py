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

  def test_read_smiles_via_context_reader(self):
    tmpdir = tempfile.mkdtemp(dir=absltest.TEST_TMPDIR.value)
    fname = os.path.join(tmpdir, "input.smi")
    with open(fname, "w") as writer:
      writer.write(SMILES)

    molecules_read = 0
    with ReaderContext(fname, FileType.SMI) as reader:
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

    set_display_smiles_interpretation_error_messages(0)
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
    stem = os.path.join(tmpdir, "somefile")
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

    with ContextWriter(stem, FileType.SMI) as writer:
      for m in mols:
        writer.write(m)

    fname = os.path.join(tmpdir, "somefile.smi")
    molecules_read = 0
    with ReaderContext(fname) as reader:
      for mol in reader:
        molecules_read += 1

    self.assertEqual(molecules_read, 10)


if __name__ == '__main__':
  #app.run(absltest.main)
  absltest.main()
