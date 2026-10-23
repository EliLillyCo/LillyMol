# Substituent Identification

The `SubstituentIdentificationLookup` binding exposes the lookup/generation side
of the command-line [`substituent_identification`](/docs/Molecule_Tools/substituent_identification.md)
tool. The database is still built by the command-line tool; Python opens that
BerkeleyDB database and asks for replacement molecules for an input `Molecule`.

This is useful when substituent replacement should be part of a Python workflow,
where generated products are inspected, filtered, scored, or passed to another
model without writing an intermediate smiles stream.

## Build a Database

Build the database exactly as for the command-line workflow:

```shell
substituent_identification -d chembl.sbs.bdb -B -R 3 -w 10 -M 10 -v \
    -Y dbproto -Y rpt=10000 chembl.smi
```

The database stores substituent fragments keyed by the EC fingerprint of the
context from which the fragment was removed. The `-R`, `-w`, `-m`, and `-M`
settings determine what contexts and fragments are available later. See the
command-line documentation for the full interpretation of these options.

## Python Lookup

```python
from lillymol import MolFromSmiles
from lillymol_bdb import SubstituentIdentificationLookup

mol = MolFromSmiles("C methane")
if mol is None:
    raise ValueError("Invalid SMILES")

with SubstituentIdentificationLookup() as lookup:
    if not lookup.open_database("/path/to/chembl.sbs.bdb"):
        raise ValueError("Cannot open substituent database")

    # Consider atoms with implicit hydrogens as growth points, like the CLI -y
    # style default-starting-point mode.
    lookup.set_default_new_molecule_starting_points(True)

    # Optional filters corresponding to common command-line controls.
    lookup.set_min_shell_radius(2)
    lookup.set_min_substituent_size(1)
    lookup.set_max_substituent_size(10)
    lookup.set_max_atoms_in_product(45)
    lookup.set_min_examples_needed_for_addition(3)
    lookup.set_max_matches_per_input_molecule(1000)

    products = lookup.generate_replacements(mol)

for product in products:
    print(product.smiles, product.name, product.donor, product.radius, product.examples)
```

Use a context manager, or call `lookup.close()`, so BerkeleyDB file handles are
released at a predictable time. This matters most when using temporary databases
or filesystems where deleting an open database file is problematic.

## Product Data

`generate_replacements(mol)` returns a list of `SubstituentReplacement` objects.
Each object contains:

| Field | Meaning |
| ----- | ------- |
| `molecule` | generated `Molecule` object |
| `smiles` | smiles for the generated product |
| `name` | name of the input molecule |
| `donor` | identifier of the molecule that supplied the replacement fragment |
| `radius` | shell radius at which the context matched |
| `examples` | number of database examples supporting the replacement fragment |
| `fragment_lost` | fragment removed from the input molecule |
| `fragment_added` | replacement fragment from the database |

The returned `Molecule` objects are owned Python objects. Mutating them does not
change the input molecule or the lookup database.

## Choosing Starting Points

There are two common lookup modes.

For growth from atoms with hydrogens, enable default starting points:

```python
lookup.set_default_new_molecule_starting_points(True)
products = lookup.generate_replacements(mol)
```

For query-driven replacement, add one or more SMARTS queries. If
`set_break_molecule_at_first_two_matched_atoms(True)` is enabled, the bond
between the first two matched query atoms is broken; the fragment containing the
second matched atom is removed and replaced.

```python
lookup.add_query_from_smarts("cO[CD1]")
lookup.set_break_molecule_at_first_two_matched_atoms(True)
products = lookup.generate_replacements(mol)
```

The query object is owned by the lookup object, so query settings are not shared
with other substructure searches.

## Notes

`add_database(dbname)` is also available as an alias for `open_database(dbname)`
when multiple databases should be searched.

The Python API intentionally covers the lookup/generation phase first. Database
construction remains a command-line operation because it is usually a large,
batch-oriented job and the command-line tool already provides progress reporting
and build-specific controls.

Some setter names mirror internal command-line concepts. The commonly useful
ones are `set_min_shell_radius`, `set_min_substituent_size`,
`set_max_substituent_size`, `set_max_atoms_in_product`,
`set_min_examples_needed_for_addition`, and
`set_max_matches_per_input_molecule`, and its alias
`set_max_molecules_per_input_molecule`.
