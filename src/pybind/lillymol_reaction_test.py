# Tests for LillyMol Reactions

from absl import app
from absl import logging
from absl.testing import absltest

from lillymol import *
from lillymol_query import *
from lillymol_reaction import *

class TestLillyMol(absltest.TestCase):
  def test_rdkit_cookbook(self):
    core = MolFromSmiles("*c1c(C)cccc1O")
    sidechain = MolFromSmiles("CN*")

    set_smirks_lost_atom_means_remove_frgment(1)

    rxn = Reaction()
    rxn.construct_from_smirks("[c:1][#0:3].[#0:4][*:2]>>[*:1]-[*:2]")
    smc = SidechainMatchConditions();
    rxn.add_sidechain_reagent(0, sidechain, smc);

    products = rxn.perform_reaction(core, sidechain)
    self.assertEqual(len(products), 1)
    self.assertEqual(products[0].unique_smiles(), "Oc1c(NC)c(C)ccc1")

if __name__ == '__main__':
  absltest.main()
