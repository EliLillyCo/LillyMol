# Hybridisation

LillyMol has a small, computed-on-demand hybridisation helper. It was added
primarily for Python code being translated from RDKit, where atom hybridisation
is often queried as a molecular property.

This functionality is intentionally modest. It is not stored on `Molecule` or
`Atom`, and it is not available in substructure searching. It should be treated
as a convenience descriptor for atoms in an instantiated molecule, not as a new
chemical perception layer.

## Python Interface

The Python binding exposes both a `Molecule` method and a free function.

```python
from lillymol import MolFromSmiles, Hybridization, hybridization, hybridization_name

mol = MolFromSmiles("CC(=O)N")

print(mol.hybridization(3))                         # Hybridization.SP2
print(hybridization(mol, 3) == Hybridization.SP2)   # True
print(hybridization_name(mol.hybridization(3)))     # SP2
```

The enum values are:

| Value | Meaning |
| ----- | ------- |
| `Hybridization.UNSPECIFIED` | no assignment |
| `Hybridization.S` | hydrogen-like `s` state |
| `Hybridization.SP` | `sp` |
| `Hybridization.SP2` | `sp2` |
| `Hybridization.SP3` | `sp3` |
| `Hybridization.SP2D` | `sp2d` |
| `Hybridization.SP3D` | `sp3d` |
| `Hybridization.SP3D2` | `sp3d2` |
| `Hybridization.OTHER` | recognised but not otherwise classified |

Invalid atom numbers raise `ValueError` in Python.

## Current Rules

The implementation first handles simple structural cases and then applies a
small number of RDKit-like exceptions. In broad terms:

* invalid atoms and atomic number zero are `UNSPECIFIED`;
* hydrogen is `S`;
* aromatic atoms are `SP2`;
* positively charged four-connected nitrogen is `SP3`;
* ordinary triple-bonded atoms are `SP`;
* ordinary double-bonded atoms are `SP2`;
* otherwise one- through four-connected atoms are `SP3`;
* five-connected atoms are `SP3D`;
* six-connected atoms are `SP3D2`.

The RDKit-alignment cases currently include:

* amide nitrogen is `SP2`;
* ester oxygen is `SP2`;
* phenolic oxygen is `SP2`;
* oxygen in `O-N=*` environments is `SP2`;
* sulfonamide nitrogen is usually `SP3`;
* a two-connected sulfonamide nitrogen attached to an aromatic carbon is `SP2`;
* four-connected phosphorus and sulfur are `SP3`, even when represented with
  `P=O` or `S=O` bonds;
* three-connected sulfur with at least one unsaturated bond is `SP3`.

These rules were chosen from comparison against RDKit on a sample of ChEMBL
molecules. They are pragmatic compatibility rules, not a promise that LillyMol
and RDKit will always assign the same hybridisation.

## Known Limits

Hybridisation perception around nitrogen can be subtle. Ring nitrogens, charged
nitrogens, sulfonamides, amides, anilines, and mixed environments can interact in
ways where RDKit makes specialised assignments. LillyMol handles the common cases
above, but deliberately avoids trying to reproduce every RDKit nitrogen edge
case.

Aromaticity definitions introduce another level of difference between implementations.

## C++ Interface

The C++ entry point is declared in `Molecule_Lib/hybridization.h`.

```c++
Hybridization HybridizationState(Molecule& mol, atom_number_t atom);
const char* ToString(Hybridization hybridization);
```

The helper computes the value when called. It does not add state to `Molecule`.
