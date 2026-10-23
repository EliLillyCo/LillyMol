#!/usr/bin/env python

import logging
import sys

from lillymol import *
from lillymol_io import *
from lillymol_query import *
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
from rdkit.Chem.Scaffolds import MurckoScaffold


# Check LillyMol and RDKit functionality

def do_lillymol_rdkit_compare(mols_lilly):
  mols_rdkit = [Chem.MolFromSmiles(m.smiles(), sanitize=True) for m in mols_lilly]


  benzene_llymol = QueryFromSmarts("c1ccccc1 benzene")
  pyridine_llymol = QueryFromSmarts("c1cnccc1 pyridine")
  # RDKit does this by default.
  benzene_llymol.set_find_unique_embeddings_only(True)
  pyridine_llymol.set_find_unique_embeddings_only(True)

  benzene_rdkit = Chem.MolFromSmarts("c1ccccc1")
  pyridine_rdkit = Chem.MolFromSmarts("c1cnccc1")

  aromaticity_differences = 0

  for l,r in zip(mols_lilly, mols_rdkit):
    if l.natoms() != r.GetNumAtoms():
      logging.error("Atom count mismatch %s", l.name())

    if abs(l.amw() - Descriptors.MolWt(r)) > 0.1:
      logging.error("AMW mismatch %f vs %f in %s", l.amw(), Descriptors.MolWt(r), l.name())

    if abs(l.exact_mass() - Descriptors.ExactMolWt(r)) > 0.001:
      logging.error("ExactMass mismatch %f vs %f in %s", l.exact_mass(), Descriptors.ExactMolWt(r), l.name())

#   Rdkit uses symmetric SSSR so comparisons not possible.
#   if l.nrings() != len(rdmolops.GetSymmSSSR(r)):
#     logging.error("Nrings mismatch %d vs %d in %s", l.nrings(), len(rdmolops.GetSymmSSSR(r)), l.name())
    
    for i in range(l.natoms()):
      if l.is_ring_atom(i) and not r.GetAtomWithIdx(i).IsInRing():
        logging.error("Ring atom mismatch atom %d in %s", i, l.name())

#   for i in range(l.natoms());
#     lsize = l.
#     if l. and not r.GetRingInfo().minAtomRingSize()

    rsize = 5
    for i in range(l.natoms()):
      if l.in_ring_of_size(i, rsize) and not r.GetAtomWithIdx(i).IsInRingSize(rsize):
        logging.warning("Ring size mismatch size %d in %s", rsize, l.name())

    # Explicit Hydrogen addition
    initial_r_usmi = Chem.MolToSmiles(r, canonical=True)
    initial_l_usmi = l.unique_smiles()

    r_with_hydrogens = Chem.AddHs(r)
    l.make_implicit_hydrogens_explicit()   # l.AddHs() does same thing.
    if l.natoms() != r_with_hydrogens.GetNumAtoms():
      logging.warning("Explicit Hydrogen count mismatch %s %d %d", l.name(),
                l.natoms(), r_with_hydrogens.GetNumAtoms())

    l.remove_all(1)   # Remove all atoms with atomic number 1.
    rback = Chem.RemoveHs(r_with_hydrogens)

    if Chem.MolToSmiles(rback, canonical=True) != initial_r_usmi:
      logging.error("Unstable rdkit H addition/removal usmi %s %s %s", l.name(),
                initial_r_usmi, Chem.MolToSmiles(rback, canonical=True))

    if l.unique_smiles() != initial_l_usmi:
      logging.error("Unstable LillyMol H addition/removal usmi %s", l.name())

    # Aromaticity differences are to be expected.
    aromaticity_diff_this_molecule = False
    for i in range(l.natoms()):
      if l.is_aromatic(i) and not r.GetAtomWithIdx(i).GetIsAromatic():
        logging.info("Aromaticity difference in %s %s %r %r", l.name(), l.smarts_equivalent_for_atom(i),
                  l.is_aromatic(i), r.GetAtomWithIdx(i).GetIsAromatic())
        aromaticity_diff_this_molecule = True
    if aromaticity_diff_this_molecule:
      aromaticity_differences += 1

    if benzene_llymol in l and not r.HasSubstructMatch(benzene_rdkit):
      logging.error("Benzene substructure match mismatch %s %r %r", l.name(),
                     benzene_llymol in l, r.HasSubstructMatch(benzene_rdkit))

    if benzene_llymol.substructure_search(l) != len(r.GetSubstructMatches(benzene_rdkit)):
      logging.warning("Benzene count mismatch in %s %d vs %d", l.name(),
                benzene_llymol.substructure_search(l),
                len(r.GetSubstructMatches(benzene_rdkit)))

    # Aromaticity differences are to be expected.
    if pyridine_llymol in l and not r.HasSubstructMatch(pyridine_rdkit):
      logging.error("Pyridine substructure match mismatch %s %r %r", l.name(),
                     pyridine_llymol in l, r.HasSubstructMatch(pyridine_rdkit))

    if pyridine_llymol.substructure_search(l) != len(r.GetSubstructMatches(pyridine_rdkit)):
      logging.warning("Pyridine count mismatch in %s %d vs %d", l.name(),
                pyridine_llymol.substructure_search(l),
                len(r.GetSubstructMatches(pyridine_rdkit)))

    scaffold_l = l.scaffold()
    scaffold_r = MurckoScaffold.GetScaffoldForMol(r)
    if scaffold_l.natoms() != scaffold_r.GetNumAtoms():
      logging.info("Atoms in scaffold mismatch %s %d vs %d",
                   scaffold_l.natoms(), scaffold_r.GetNumAtoms())


  logging.info("Found aromaticity differences in %d molecules", aromaticity_differences)

def lillymol_rdkit_compare(argv):
  if len(argv) == 1:
    logging.error("Must specify input file to use")

  mols = slurp(argv[1])
  print(f"Read {len(mols)} molecules from {argv[1]}", file=sys.stderr)

  do_lillymol_rdkit_compare(mols)

if __name__ == '__main__':
  lillymol_rdkit_compare(sys.argv)
