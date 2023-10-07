from absl import app
from absl import flags

"""Demo app to show LillyMol python bindings.
  Convert 4-pyidol to 4-pyidone forms.
  This was done in order to explore a new chemical standardisation.
  This implementation was done in python, and when it was found to
  be doing what we wanted, it was translated to C++.

  Run time for processing the corporate collection with this file
  was about 6.5 minutes. Doing the same thing with the c++ version
  was about 1 minute 15 seconds.

  For some tasks this trade-off may be attractive.
"""

from lillymol import *

def change_to_pyridone(m:Molecule, ring:Ring, n_index:int, oh_index:int, oh:int)->None:
  """Convert `ring` to 4-pyridone form
    Args:
      m: molecule
      r: 6 membered aromatic ring
      n_index: index in r of [nD2]
      oh_index: index in r of c-[OH]
      oh: atom number of [OH]

    Note that this only works because python has a copy of the rings.
    Once a bond is changed, rings are destroyed. But safe here because
    the ring is not changing.
  """

  # The para exocyclic -O becomes =O
  m.set_bond_type_between_atoms(ring[oh_index], oh, BondType.DOUBLE_BOND)

  # Introduce sequence of alternating single and double bonds into `ring`,
  # beginning at the Nitrogen.
  to_place = [1, 2, 1, 1, 2, 1]
  prev = n_index
  for i,ndx in enumerate(range(n_index + 1, n_index + 7)):
    if ndx > 5:
      ndx = ndx % 6
    if to_place[i] == 1:
      m.set_bond_type_between_atoms(ring[prev], ring[ndx], BondType.SINGLE_BOND)
    else:
      m.set_bond_type_between_atoms(ring[prev], ring[ndx], BondType.DOUBLE_BOND)
    prev = ndx

  print(m.smiles() + ' ' + m.name())
    

  return

def ring_is_four_pyridone(m:Molecule, ring:Ring)->None:
  """ See if `ring` is a 4-pyridone"""
  # The indices of things in the ring, or the OH
  n_index = -1
  oh_index = -1
  oh = -1;

  # Look for [n] and OH in the ring
  for ndx,a in enumerate(ring):
    if m.atomic_number(a) == 7:
      if n_index >= 0:
        return
      n_index = ndx
      continue

    # No other heteroatoms
    if m.atomic_number(a) != 6:
      continue

    if m.ncon(a) == 2:
      continue

    # Look for exocyclic double bonds and the [OD1]
    for bond in m[a]:
      o = bond.other(a);
      if o in ring:
        continue

      if not bond.is_single_bond():
        return

      if m.atomic_number(o) != 8:
        continue
      if m.ncon(o) != 1:
        continue
      if oh_index >= 0:
        return
      oh_index = ndx
      oh = o

  if n_index < 0 or oh_index < 0:
    return

  # print(f'{ring} N {n_index} c-OH {oh_index} OH {oh}')
  # The two indices must be 3 apart in `ring`. 
  if (n_index + 3) % 6 == oh_index:
    pass
  elif (oh_index + 3) == n_index:
    pass
  else:
    return

  change_to_pyridone(m, ring, n_index, oh_index, oh)

def four_pyridone(m:Molecule):
  m.compute_aromaticity_if_needed();
  for r in m.rings():
    if r.size() != 6:
      continue
    if not r.is_aromatic():
      continue

    ring_is_four_pyridone(m, r)

def main(argv):
  """Chemical standardisation expt"""
  with open(argv[1], 'r') as reader:
    for line in reader:
      m = MolFromSmiles(line)
      if m is None:
        print(f'Skipping {line.rstrip()}', file=sys.stderr)
        continue
      four_pyridone(m)

if __name__ == '__main__':
  app.run(main)
