import sys
from itertools import combinations
import logging

from rdkit import Chem

def attach_atom(mol_in, atomic_symbol, smarts, num_sub, seen):
    pt = Chem.GetPeriodicTable()
    atomic_num = pt.GetAtomicNumber(atomic_symbol)
    out_mol_list = []
    query = Chem.MolFromSmarts(smarts)
    sresult = mol_in.GetSubstructMatches(query)
    if sresult is None:
      return out_mol_list
    match_atms = [x[0] for x in sresult]
    n_combos = combinations(match_atms, num_sub)
    for combo in n_combos:
        new_mol = Chem.RWMol(mol_in)
        for idx in combo:
            new_idx = new_mol.AddAtom(Chem.Atom(atomic_num))
            new_mol.AddBond(idx, new_idx, order=Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(new_mol)
        smi = Chem.MolToSmiles(new_mol)
        if smi not in seen:
            seen.add(smi)
            out_mol_list.append(new_mol)
    return out_mol_list

def main(argv):

  seen = set()
  atomic_symbol = "F"
  smarts = "[cH]"
  num_sub = 2

  with open(argv[1], 'r') as reader:
    for line in reader:
      mol = Chem.MolFromSmiles(line)
      if mol is None:
        continue
      result = attach_atom(mol, atomic_symbol, smarts, num_sub, seen)
      seq = 0
      for x in result:
        print(f"{Chem.MolToSmiles(x)} #{x.GetProp('_Name')}.{seq}")
        seq += 1


if __name__ == '__main__':
  main(sys.argv)

