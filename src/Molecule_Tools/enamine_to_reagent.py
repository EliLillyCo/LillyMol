#!/usr/bin/env python

# Input file looks like
# [U]CCC        id      1,2,3   rxn     notused
# Convert to a recordio file of ReagentData serialized protos.

from absl import app
from absl import logging

from lillymol import *

from Molecule_Tools.reagent_substructure_search_pb2 import *

def transuranium_to_isotopes(mol:Molecule)->int:
  """Remove U, Np, Pu atoms in `mol` and place isotopes on the
     attached atoms. 
     Isotope varies with the kind of bond removed.
     Return the number of atoms removed.
  """
  rc = 0
  for i in range(mol.natoms() - 1, 0, -1):
    atom = mol[i]
    if atom.atomic_number() < 92:
      continue

    bond = atom[0];
    connected = bond.other(i);
    if bond.is_single_bond():
      mol.set_isotope(connected, 1)
    elif bond.is_double_bond():
      mol.set_isotope(connected, 2)
    elif bond.is_triple_bond():
      mol.set_isotope(connected, 3)

    mol.remove_atom(i)
    rc += 1

  return rc;

def enamine_to_reagent_2(line, usmi2reagent):
  f = line.split()

  mol = MolFromSmiles(f[0])
  if mol is None:
    return

  ncon = transuranium_to_isotopes(mol)

  if mol.unique_smiles() in usmi2reagent:
    proto = usmi2reagent.get(mol.unique_smiles())
  else:
    proto = ReagentData()
    proto.smiles = mol.unique_smiles()
    proto.ncon = ncon
    proto.natoms = mol.natoms()
    usmi2reagent[mol.unique_smiles()] = proto

  new_reaction = proto.reaction.add()
  new_reaction.name = f[3]
  new_reaction.x = int(f[2])

def enamine_to_reagent(argv):
  if len(argv) == 1:
    logging.error("Must specify enamine reagent file")
    return 1

  output_file = 'enamine.recordio'

  # A mapping from unique smiles to a ReagentData
  usmi2reagent = {}

  with open(argv[1], "r") as input:
    for ndx, line in enumerate(input):
      enamine_to_reagent_2(line, usmi2reagent)
      if ndx % 10000 == 0:
        logging.info("Processed %d records", ndx)

  logging.info("Identified %u reagents", len(usmi2reagent))

  with open(output_file, "wb") as output:
    for proto in usmi2reagent.values():
      s = proto.SerializeToString()
      slen = f"{len(s)}\n"
      output.write(slen.encode('utf-8'))
      output.write(s)
      

if __name__ == '__main__':
  app.run(enamine_to_reagent)
