#!/usr/bin/env python3

# parse the KipKake output from numbered_smiles
# O([1C]1([30H])[2C]([31H])([32H])[3N]([4C]1([33H])[34H])[5S](=[6O])(=[7O])[8C]([35H])([36H])[9C]([37H])([38H])[10O][11C]([39H])([40H])[41H])[12C]1=[13C]([42H])[14C](=[15C]([16C](=[17C]1[45H])[44H])[18N]([46H])[19C](=[20O])[21C@]1([47H])[22C@]([48H])([23C]1([49H])[50H])[24C]1=[25C]([51H])[26C](=[27C]([53H])[28N]=[29C]1[54H])[52H])[43H]
# 3186120 idx:28 pka_type:basic pka:4.812224

from absl import app
from absl import flags

from lillymol import *

def from_qupkake_line(line):
  mol = MolFromSmiles(line)
  name = mol.name()
  f = name.split(' ')
  id = f[0]
  idx = f[1].split(':')[1]
  pka_type = f[2].split(':')[1]
  pka = f[3].split(':')[1]
  idx = int(idx)
  pka = float(pka)

  x = mol.first_atom_with_isotope(idx)

  mol.remove_isotopes()

#  print(f"ELE {mol.atomic_number(idx)}")
  mol.set_isotope(x, int(1000.0*pka))

  mol.remove_all(1)
  print(f"{mol.smiles()} {id} {pka_type}")


def from_qupkake(argv):
  with open(argv[1], 'r') as file:
    for line in file:
      from_qupkake_line(line)
  

if __name__ == '__main__':
  app.run(from_qupkake)
