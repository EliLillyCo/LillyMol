from absl import app
from absl import logging

from lillymol import *
from lillymol_io import *
from lillymol_query import *
from lillymol_reaction import *

# Demonstrations of reactions in LillyMol python.
# Some reagent sets to be used. These could also be in files.

ACIDS="""OC(=O)C acid1
OC(=O)CC acid2
OC(=O)CCC acid3"""

AMINES="""NC amine1
NCC amine2
NCCC amine3"""

PROTECTED_AMINES="""CC(C)(C)OC(=O)NC tb1,
CC(C)(C)OC(=O)NCC tb2
CC(C)(C)OC(=O)NCCC tb3
c12ccccc1c3ccccc3C2COC(=O)NC fmc1
c12ccccc1c3ccccc3C2COC(=O)NCC fmc2
c12ccccc1c3ccccc3C2COC(=O)NCCC fmc3"""

# Reactions that process the reagents above

# Both tboc and fmoc are removed.
# Note the use of multiple smarts.
# The important factor is that in both smarts
# the bond to be broken must be the same indices.
DEPROTECT_RXN="""name: "deprotection"
scaffold {
  smarts: "NC(=O)-OC(C)(C)C"
  smarts: "NC(=O)-OCC1c2ccccc2c2ccccc21"
  break_bond {
    a1: 0
    a2: 1
  }
  remove_fragment: 1
  isotope {
    atom: 0
    isotope: 1
  }
  match_conditions {
    find_unique_embeddings: true
  }
}"""

# Add t-Boc to am amine
ADD_FIXED_REAGENT_RXN="""name: "t-Boc protect amine"
scaffold {
  id: 0
  smarts: "[ND1]-[CX4]"
}
sidechain {
  id: 1
  reagent: "C(=O)OC(C)(C)C"
  smarts: "C(=O)OC(C)(C)C"
  join {
    a1: 0
    a2: 0
  }
}"""


ACID_AMINE="""name: "acid + amine"
scaffold {
  smarts: "[OD1]-[CD3]=O"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "[NH2]-[CX4]"
  join {
    a1: 1
    a2: 0
  }
}
"""

def get_molecules(many_smiles) ->list[Molecule]:
  """many_smiles is a newline separated list of smiles.
     return a list of the corresponding molecules.
     Errors are skipped
  """
  result = []
  for smiles in many_smiles.split('\n'):
    mol = MolFromSmiles(smiles)
    if mol is None:
      logging.error("Invalid smiles %s", smiles)
  
    result.append(mol)

  return result

def form_name(scaffold, rxn, iter, product):
  """Assign the name of `product` given that it has been generated
    by `rxn` from `scaffold` and with `iter` olding information about
    sidechains used
  """
  names = rxn.reagent_names(iter)
  new_name = scaffold.name()
  for i in range(1, len(names)):
    new_name += f" + {names[i]}"
  product.set_name(new_name)
  return

def add_fixed_reagent():
  """A reaction with a single sidechain.
    A t-Boc group is added to some amines
  """
  amines = get_molecules(AMINES)

  rxn = Reaction()
  if not rxn.construct_from_textproto(ADD_FIXED_REAGENT_RXN):
    logging.error("Invalid reaction %s", ADD_FIXED_REAGENT_RXN)

  for amine in amines:
    if rxn.in_place_transformations(amine):
      logging.info("%s protected amine", amine.smiles())
    else:
      logging.error("%s no t-Boc protection", amine.smiles())

def in_place_transformation():
  """Make changes to a molecule in place, no sidechains added.
  """
  protected_amines = get_molecules(PROTECTED_AMINES)

  rxn = Reaction()
  if not rxn.construct_from_textproto(DEPROTECT_RXN):
    logging.error("Invalid reaction %s", DEPROTECT_RXN)

  for amine in protected_amines:
    if not rxn.in_place_transformations(amine):
      logging.error("%d not deprotected", amine.name())
    logging.info("%s deprotected amine", amine.smiles())

def full_enumeration():
  # Enumerate a trxn-like acid_amine reaction

  rxn = Reaction()
  if not rxn.construct_from_textproto(ACID_AMINE):
    logging.error("Invalid reaction %s", ACID_AMINE)

  logging.info("Reaction %s has %d sidechains", rxn.name(), rxn.number_sidechains())

  scaffolds = get_molecules(ACIDS)
  sidechains = get_molecules(AMINES)

  logging.info("Read %d scaffolds and %d sidechains", len(scaffolds), rxn.number_sidechains())

  # Needed if there are multiple substructure matches in a sidechain.
  sidechain_match_conditions = SidechainMatchConditions()

  for sidechain in sidechains:
    if not rxn.add_sidechain_reagent(0, sidechain, sidechain_match_conditions):
      logging.error("Invalid sidechain %s, adjust settings in sidechain_match_conditions", smiles)

  # As scaffolds are processed, we iterate through sidechains.
  iter = ReactionIterator(rxn)

  for scaffold in scaffolds:
    matches = rxn.substructure_search_matches(scaffold)
    if not matches:
      logging.error("No scaffold match to %s matched %d atoms", scaffold.name(), rxn.max_query_atoms_matched_in_search())
      return

    iter.reset()
    while iter.active():
      for match in matches:
        product = rxn.perform_reaction(scaffold, match, iter)
        if product is None:
          logging.error("No reaction from %s with %s", scaffold.name(), rxn.name())
          continue

        form_name(scaffold, rxn, iter, product)
        logging.info("made %s %s", product.smiles(), product.name())
        iter.increment()

def reaction_examples(argv):
  in_place_transformation()
  add_fixed_reagent()
  full_enumeration()

if __name__ == "__main__":
  app.run(reaction_examples)
