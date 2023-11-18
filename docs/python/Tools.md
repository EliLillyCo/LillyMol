# Tools

The python bindings come with a small number of tools, mostly derived
from tools otherwise avaialble as c++ executables.

## Position3D
This is also covered in the [trxn](/docs/Molecule_Tools/trxn.md) documentation.

A common operation in 3d is to react an arbitrary fragment with a core that
is fixed. Exploring different substituents at given position(s) on a scaffold.
In order to do this, the fragment must be translated and rotated so
that the atoms across which the new bond will be formed are facing each other.

It is necessary to identify the atom in the scaffold and the atom in the
fragment that will form the bond. There are `a1` and `a2` in the pseudo code
below.

The important thing to note is that `a2` must be converted from an atom number
in `fragment` to an atom number in the merged molecule, so it is added
to the number of atoms in the scaffold - which could also have been
derived from `scaffold.natoms()`.

We make a copy of `scaffold` because it will be changed and we want to
re-use it across other fragments.
```
a1 = 4
a2 = len(scaffold) + 5

mcopy = copy.copy(scaffold)
mcopy += fragment
Position3D(mcopy, a1, 1.35, a2)
mcopy.add_bond(a1, a2, BondType.SINGLE_BOND)
```
At this stage, we should have a 3d plausible arrangement of atoms `a1` and
`a2`, separated by 1.35 Angstroms.

The next task might be to perform a conformational search around the newly
formed bond. This could of course be done by successively setting the
dihedral angle, but that would be inefficient. Instead we can do a one-time
generation of all coordinates implied by a conformational scan. 
```
angle = 30
bump_check = 2.5
coords = mcopy.dihedral_scan(a1, a2, angle, bump_check)

for conf in coords:
  mcopy.set_coordindates(conf)
```
which privides a very cheap means of generating conformers. The number of
conformers returned will be influenced by the value assigned to `bump_check`
which specifies a minimum distance separating any moving or non moving atom.
Now, of course these conformers are entirely naieve, they know nothing about
their energetic implications beyond passing the bump check. More expensive
methods will be needed in order to evaluate each conformer.

