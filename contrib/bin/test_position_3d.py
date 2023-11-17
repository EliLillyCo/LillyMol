# Tester for Position3D

import copy
import random

from absl import app
from absl import logging

from lillymol import *
from lillymol_io import *

def test_3d_position(mol: Molecule, writer):
  """Breaks bonds in `mol` and regenerates with Position3D.
  For each bond in `mol` that seems suitable for breaking,
  remove that bond, randomly reposition one of the resulting
  fragments, then use Position3D to re-orient the fragment
  and then re-form the bond.

  Args:
    mol: Molecle
    writer: write re-formed molecules here.
  """
  mol.remove_all_chiral_centres()
  mol.compute_aromaticity_if_needed()

  # Write the parent once if any variants are written.
  parent_written = False

  for bond in mol.bonds():
    if bond.nrings() > 0:
      continue

    a1 = bond.a1()
    a2 = bond.a2()
    if mol.ncon(a1) == 1 or mol.ncon(a2) == 1:
      continue

    mcopy = copy.copy(mol)

    mcopy.remove_bond_between_atoms(a1, a2)
    frag = mcopy.get_fragment_membership()
    fixed_frag = frag[a1]
    moving_frag = frag[a2]

    # Randomly shift the atoms in moving_frag
    dx = random.uniform(-10, 10)
    dy = random.uniform(-10, 10)
    dz = random.uniform(-10, 10)
    mcopy.translate(frag, moving_frag, dx, dy, dz)

    rotx = random.uniform(-1, 1)
    roty = random.uniform(-1, 1)
    rotz = random.uniform(-1, 1)
    axis = Coordinates(rotx, roty, rotz)
    angle = random.uniform(-3.14, 3.14)
    axis.normalise()
    mcopy.rotate(frag, moving_frag, axis, angle)

    #mcopy.set_name("rotated")
    #writer.write(mcopy)

    # print(f"Before positioning {mcopy.distance_between_atoms(a1, a2)}")
    Position3D(mcopy, a1, 2.0, a2)
    # print(f"After positioning {mcopy.distance_between_atoms(a1, a2)}")
    mcopy.add_bond(a1, a2, BondType.SINGLE_BOND)
    if not parent_written:
      writer.write(mol)
      parent_written = True

    mcopy.set_name("translated")
    writer.write(mcopy)

    angle = 45
    bump_check = 3.0
    coords = mcopy.dihedral_scan(a1, a2, angle, bump_check)
    logging.info("%s generated %d conformations", mol.name(), len(coords))

    for conf in coords:
      mcopy.set_coordinates(conf)

def main(argv):
  """Test Position3D
  """
  if len(argv) == 1:
    logging.error("Must specify input file")

  mols = slurp(argv[1])
  logging.info("Contains %d molecules", len(mols))
  with ContextWriter("reformed.sdf", FileType.SDF) as writer:
    with ReaderContext(argv[1]) as reader:
      for mol in reader:
        test_3d_position(mol, writer)

if __name__ == '__main__':
  app.run(main)
