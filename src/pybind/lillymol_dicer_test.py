# Tests for dicer
import os
import tempfile

from absl import app
from absl import logging
from absl.testing import absltest

from google.protobuf import text_format

from lillymol import *
from lillymol_tools import *

class TestDicer(absltest.TestCase):
  def test_c4_no_break_cc(self):
    dicer = Dicer()
    m = MolFromSmiles("CCCC")
    frags = dicer.dice(m)
    self.assertEqual(len(frags), 0)

  def test_c4_break_cc(self):
    dicer = Dicer()
    dicer.set_break_cc_bonds(True)
    m = MolFromSmiles("CCCC")
    frags = dicer.dice(m)
    self.assertDictEqual(frags, {'CC': 2, 'CCC': 2, 'C': 2})

  def test_1313430(self):
    dicer = Dicer()
    m = MolFromSmiles("NOC1CCC(N)CC1 CHEMBL1213430")
    frags = dicer.dice(m)
    self.assertDictEqual(frags, {'N': 2, 'NC1CCCCC1': 1, 'NOC1CCCCC1': 1, 'OC1CCC(N)CC1': 1, 'ON': 1})

  def test_1313430_iso(self):
    dicer = Dicer()
    dicer.set_label_join_points(1)
    m = MolFromSmiles("NOC1CCC(N)CC1 CHEMBL1213430")
    frags = dicer.dice(m)
    self.assertDictEqual(frags, {'NOC1CC[1CH2]CC1': 1, 'NC1CC[1CH2]CC1': 1, '[1OH]N': 1, '[1OH]C1CCC(N)CC1': 1, '[1NH3]': 2})

  def test_2354634Recap(self):
    dicer = Dicer()
    dicer.set_label_join_points(1)
    dicer.set_work_like_recap(True)
    m = MolFromSmiles("O=C(N(C)C1=C(N=C2N1C=C(C=C2)C(=O)NCCOC1=CC=C(OC)C=C1)CC)CC1=CC=CC=C1 CHEMBL2354634")
    frags = dicer.dice(m)
    self.assertDictEqual(frags, {'O=[1CH][1NH2]': 2, '[1CH3]C': 1, '[1CH3][1CH3]': 1, '[1CH4]': 3, '[1OH2]': 2, '[1cH]1cc[1cH]cc1': 1, '[1cH]1ccccc1': 1, '[n]1c2[n](c[1cH]cc2)[1cH][1cH]1': 1})

if __name__ == '__main__':
  #app.run(absltest.main)
  absltest.main()
