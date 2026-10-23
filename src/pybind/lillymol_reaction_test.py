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

  def test_simple_multiple_reagents(self):
    reagents = []
    reagents.append(MolFromSmiles("O-C(=O)c1ccc(Cl)cc1 scaffold"))
    reagents.append(MolFromSmiles("Nc1ccc(S)cc1 R1"))
    reagents.append(MolFromSmiles("C R2"))
    rxn_textproto = """scaffold {
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
    rxn = Reaction();
    self.assertTrue(rxn.construct_from_textproto(rxn_textproto))
    product = rxn.perform_reaction(reagents)
    self.assertIsNotNone(product)
    self.assertEqual(product.unique_smiles(), "O=C(Nc1ccc(C)cc1)c1ccccc1")

if __name__ == '__main__':
  absltest.main()
