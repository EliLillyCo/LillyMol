# 2D depiction

The optional `lillymol_depict` module generates two-dimensional coordinates for
a LillyMol `Molecule`. It uses Schrodinger's coordgen library and is built when
`BUILD_SCHRODINGER_2D` is enabled. The main `lillymol` module does not depend on
coordgen.

## Generate coordinates in place

```python
import lillymol
import lillymol_depict as depict

mol = lillymol.MolFromSmiles("CCOc1ccccc1")
result = depict.generate_2d_coordinates(mol)
if result is None:
    raise RuntimeError("Could not generate a usable layout")

print(result.closest_nonbonded_approach)
```

The function replaces the molecule's coordinates, sets every z coordinate to
zero, and returns a `Coords2DResult`. It returns `None` on failure and leaves the
original coordinates untouched.

Options can be passed directly as keyword arguments:

```python
result = depict.generate_2d_coordinates(
    mol,
    bond_length=1.5,
    precision=depict.Precision.BEST,
    skip_minimization=False,
    even_angles=False,
    centre=True,
    honour_cis_trans=True,
)
```

For repeated calls, construct an options object:

```python
options = depict.Coords2DOptions(
    precision=depict.Precision.STANDARD,
    honour_cis_trans=True,
)

for mol in molecules:
    result = depict.generate_2d_coordinates(mol, options)
```

`Precision` has `QUICK`, `STANDARD`, and `BEST` values. Higher precision spends
more time resolving difficult layouts; most ordinary molecules produce the
same result at all three settings.

The result reports:

- `closest_nonbonded_approach`: the closest non-bonded atom distance, or a
  negative value when the molecule has fewer than three atoms;
- `cis_trans_bonds`: configured double bonds that constrained the layout;
- `cis_trans_honoured`: constrained bonds drawn with the requested geometry.

## Generate a copy

Use `generate_2d_coordinates_copy` when the input coordinates must be retained:

```python
depicted = depict.generate_2d_coordinates_copy(mol, options)
if depicted is not None:
    mol2d, result = depicted
```

It returns `(Molecule, Coords2DResult)` on success or `None` on failure. The
input molecule is never changed.

## Wedge and hash bonds

Coordinate generation and tetrahedral wedge assignment are separate operations:

```python
result = depict.generate_2d_coordinates(mol)
if result is not None:
    wedge = depict.assign_wedge_bonds(mol)
    print(wedge.chiral_centres, wedge.wedged, wedge.unresolved)
```

`assign_wedge_bonds` clears old wedge/hash annotations and assigns bonds that
reproduce the chirality already held by the molecule. Its `WedgeResult` reports
`chiral_centres`, `wedged`, `not_stereogenic`, and `unresolved` counts.

Explicit hydrogens are laid out as ordinary atoms. Remove them before coordinate
generation when the desired depiction should use implicit hydrogens.
