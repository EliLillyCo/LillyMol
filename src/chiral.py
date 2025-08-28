
from absl import app

from lillymol import *
from lillymol_io import *

def compress_halogens(z):
  if z < 35 :
    return z
  if z == 35:
    return 17

  if z == 53:
    return 17

  return 17  # huh!

def attached_score(mol, zatom, avoid):
  found = []
  for bond in mol[zatom]:
    o = bond.other(zatom)
    if o == avoid:
      continue

    z = mol.atomic_number(o)
    z = compress_halogens(z)
    if bond.is_single_bond():
      found.append(z)
    elif bond.is_double_bond():
      found.append(10 * z)
    else:
      found.append(100 * z)

  hash = 0
  for s in sorted(found):
    hash += 100 * hash + s

  return 100000 + hash
  

def main(argv):
  seen = dict()
  n4 = 0
  with ReaderContext(argv[1]) as reader:
    for mol in reader:
      for c in mol.chiral_centres():
#       print(c)

        centre = c.atom()
        attached = set()
        ok = True
        for bond in mol[centre]:
          o = bond.other(centre)
          ascore = attached_score(mol, o, centre)
          if ascore in attached:
            ok = False
          attached.add(ascore)

        if not ok:
          continue

        s = sorted(attached)
        hash = 1000000 * s[0] + 10000 * s[1] + 100 * s[2]
        if len(s) == 4:
          n4 += 1
          hash += s[3]
        if hash in seen:
          seen[hash] += 1
        else:
          seen[hash] = 1

  for hash, count in seen.items():
    print(f"{hash} - {count}")
  print(f"{n4} four connected atoms")

if __name__ == '__main__':
  app.run(main)
