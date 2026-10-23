# molecule_filter

`molecule_filter` is a fast first-pass molecule filtering tool. It reads SMILES
or any other LillyMol-supported molecule input, applies a textproto set of hard
molecular-property constraints, and writes the molecules that pass.

It can also score surviving molecules with configurable piecewise-linear
utility functions. In that mode, hard filters still run first; utility values are
computed only for molecules that pass those filters.

## Quick start

Create a textproto configuration file:

```textproto
min_natoms: 10
max_natoms: 40
min_nrings: 1
max_nrings: 6
max_alogp: 6
```

Run the filter:

```shell
molecule_filter -F config.textproto large.smi > passed.smi
```

To see a summary of why molecules were rejected:

```shell
molecule_filter -F config.textproto -v large.smi > passed.smi
```

To also keep rejected molecules:

```shell
molecule_filter -F config.textproto -B rejected large.smi > passed.smi
```

The configuration file is based on
[`molecule_filter.proto`](/src/Molecule_Tools/molecule_filter.proto).

## Input and output

### Normal filter mode

If the textproto contains only hard filters and no `utility` messages,
`molecule_filter` writes surviving molecules as normal SMILES records.

```text
smiles id
```

If the molecule is modified during processing, for example by `-l`, `-c`, or
chemical standardisation, the modified SMILES is written. Otherwise, the original
input line is written unchanged.

### Utility table mode

If the textproto contains one or more `utility` messages, hard filters are still
applied first. Molecules that fail any hard filter are discarded and no utility
values are computed for them.

By default, surviving molecules are written as a descriptor-style table:

```text
Name utility1 utility2 ... overall_utility
```

For example, this configuration:

```textproto
min_natoms: 3
utility {
  name: "natoms"
  point { x: 0 y: 0 }
  point { x: 10 y: 1 }
}
utility {
  name: "heteroatom_count"
  point { x: 0 y: 0 }
  point { x: 2 y: 1 }
}
```

applied to:

```text
CC ethane
CCO ethanol
CCCC butane
```

writes:

```text
Name natoms heteroatom_count overall_utility
ethanol 0.3 0.5 0.4
butane 0.4 0 0.2
```

`ethane` is removed by the hard `min_natoms` filter before utility values are
computed.

### Utility SMILES mode

If the textproto contains utilities and you want the output to remain a regular
SMILES file, use `-u`. In this mode, individual utility columns are suppressed
and no header is written.

```shell
molecule_filter -F config.textproto -u input.smi > scored.smi
```

Output is:

```text
smiles id overall_utility
```

For the example above:

```text
CCO ethanol 0.4
CCCC butane 0.2
```

This is useful when the result will be consumed by tools that expect a SMILES
file, while still retaining the overall utility as an extra column.

## Command-line options

| Option | Meaning |
|--------|---------|
| `-F <fname>` | Textproto configuration file. Required. |
| `-l` | Reduce each input record to the largest fragment before molecule-level processing. |
| `-c` | Remove all chirality before molecule-level processing. |
| `-g ...` | Apply LillyMol chemical standardisation options. |
| `-B <fname>` | Write rejected input records to `<fname>.smi`. |
| `-u` | With utilities, write `smiles id overall_utility` instead of a descriptor table. |
| `-v` | Verbose output, including rejection summary. |

Standard LillyMol aromaticity and element options are also accepted.

## Hard filters

Hard filters are specified directly in the `Requirements` textproto. The tool is
optimised to test cheap conditions before expensive ones. For SMILES input, atom
and ring counts can often be checked from the SMILES string before constructing a
full `Molecule` object.

A large example showing the available hard filters is:

```textproto
min_natoms: 10
max_natoms: 40
min_heteroatom_count: 1
max_heteroatom_count: 20
min_heteroatom_fraction: 0.05
max_heteroatom_fraction: 0.90
min_nrings: 1
max_nrings: 6
min_aromatic_ring_count: 1
max_aromatic_ring_count: 4
min_aliphatic_ring_count: 0
max_aliphatic_ring_count: 4
min_rotatable_bonds: 2
max_rotatable_bonds: 10
min_ring_system_size: 0
max_ring_system_size: 3
min_aromatic_rings_in_system: 0
max_aromatic_rings_in_system: 3
min_tpsa: 20
max_tpsa: 180
min_alogp: -2
max_alogp: 6
min_xlogp: -2
max_xlogp: 6
min_hba: 2
max_hba: 10
min_hbd: 0
max_hbd: 5
max_largest_ring_size: 7
exclude_non_organic: true
exclude_isotopes: true
min_halogen_count: 0
max_halogen_count: 7
max_distance: 30
min_sp3_carbon: 1
min_sp3_carbon_fraction: 0.20
max_aromatic_density: 0.8
max_chiral: 2
max_number_fragments: 2
```

As written, this example is not intended as a recommended filter. It is an
illustration of available fields. In particular, it is usually unnecessary to
filter on both ALOGP and XLOGP. XLOGP is more expensive; if speed is critical,
ALOGP may be preferred.

## Utility functions

A `utility` message names one molecular feature and defines a piecewise-linear
function that maps the raw feature value to a utility value.

```textproto
utility {
  name: "natoms"
  weight: 2
  point { x: 10 y: 0 }
  point { x: 25 y: 1 }
  point { x: 40 y: 0 }
}
```

Rules:

- At least two `point` values are required.
- Points do not need to be entered in sorted order; the tool sorts them by `x`.
- Duplicate `x` values are rejected.
- `y` values must be in `[0, 1]`.
- Values below the smallest `x` use the first `y` value.
- Values above the largest `x` use the last `y` value.
- Values between adjacent points use linear interpolation.
- If `weight` is omitted, it defaults to `1`.

### Combining utilities

When multiple utility functions are present, they are combined into
`overall_utility`. The default is weighted average.

```textproto
utility_combination: UTILITY_COMBINATION_WEIGHTED_AVERAGE
```

Available modes are:

| Mode | Meaning |
|------|---------|
| `UTILITY_COMBINATION_WEIGHTED_AVERAGE` | Default. `sum(weight * utility) / sum(weight)`. |
| `UTILITY_COMBINATION_WEIGHTED_SUM` | `sum(weight * utility)`. |
| `UTILITY_COMBINATION_PRODUCT` | Product of utility values. Weights are ignored. |
| `UTILITY_COMBINATION_MIN` | Minimum utility value. Weights are ignored. |
| `UTILITY_COMBINATION_MAX` | Maximum utility value. Weights are ignored. |

### Utility feature names

Canonical utility feature names are:

```text
natoms
nrings
heteroatom_count
heteroatom_fraction
aromatic_ring_count
aliphatic_ring_count
rotatable_bonds
max_ring_system_size
aromatic_rings_in_system
tpsa
alogp
xlogp
hba
hbd
largest_ring_size
halogen_count
max_distance
sp3_carbon
sp3_carbon_fraction
aromatic_density
chiral
number_fragments
```

Recognised aliases are:

```text
heteroatoms
aromatic_rings
aliphatic_rings
rotbond
ring_system_size
halogens
longest_path
nfrag
fragments
fraction_csp3
fcsp3
```

`sp3_carbon_fraction` is the fraction of carbon atoms that are fully saturated
sp3 carbons. The denominator is the number of carbon atoms, not the total atom
count. Molecules with no carbon atoms have `sp3_carbon_fraction` 0.0.

Feature names are validated when the configuration is read. Unknown names cause
initialisation to fail before any molecules are processed.

## Examples

### Size and ring filter

```textproto
min_natoms: 10
max_natoms: 40
min_nrings: 1
max_nrings: 6
```

```shell
molecule_filter -F config.textproto input.smi > passed.smi
```

### Largest fragment filtering

Use `-l` when salts or mixtures should be reduced before molecule-level checks.

```shell
molecule_filter -F config.textproto -l input.smi > passed.smi
```

If `max_distance` is specified and `-l` is not present, `molecule_filter`
automatically enables largest-fragment processing, because longest-path checks
are meaningful only within a fragment.

### Hard filter plus utility score

```textproto
min_natoms: 10
max_natoms: 40
max_alogp: 6
utility {
  name: "natoms"
  point { x: 10 y: 0 }
  point { x: 25 y: 1 }
  point { x: 40 y: 0 }
}
utility {
  name: "alogp"
  point { x: 1 y: 0 }
  point { x: 3 y: 1 }
  point { x: 6 y: 0 }
}
```

This first removes molecules outside the hard atom-count and ALOGP limits. Only
surviving molecules are scored with the two utility functions.

## Performance notes

The order of computation is chosen to avoid expensive work when cheaper checks
can reject a molecule first. Atom counts, ring counts, fragment counts, isotope
checks, and related simple properties are handled before more expensive
calculations such as TPSA, ALOGP, XLOGP, ...

Utility values are computed only after hard filters pass. If utilities request a
feature that was not part of the hard filters, that feature is computed lazily for
surviving molecules.

This task is usually pleasingly parallel. For very large collections,
`molecule_filter` is commonly run in sharded parallel workflows; see also
`molecule_filter_parallel.sh`.

## Rejection summary

With `-v`, `molecule_filter` writes a summary of rejection counts to stderr. For
example, processing 2.4M ChEMBL molecules might report:

```text
Read 2409270 molecules, passed 1488785 0.61794
7905 non organic
1901 isotope
6582 too few atoms 10
253993 too many atoms 40
18088 too few rings 1
10255 too many rings 6
282 too few heteroatoms 1
2170 min heteroatom fraction 0.05
13 max heteroatom fraction 0.9
41203 too few aromatic rings 1
64474 too many aromatic rings 4
0 too few aliphatic rings 0
901 too many aliphatic rings 4
22132 ring systems too large 3
0 too many aromatic rings in system 3
7964 ring too large 7
84062 too few rotatable bonds 2
60752 too many rotatable bonds 10
14284 low TPSA 20
3559 high TPSA 180
294 low ALOGP -2
56821 high ALOGP 6
13937 too few HBA 2
34722 too many HBA 10
0 too few HBD 0
13070 too many HBD 5
1761 too many halogens 7
0 molecules too long 30
625 too few CSP3 1
47640 aromatic density too high 0.8
151095 too many chiral centres 2
18800 too many fragments 2
```

Some rules may show zero matches because earlier rules already rejected the
molecules that would have failed them. The order shown in the summary follows the
proto fields; the actual computation order is implementation-defined and chosen
for speed.

## Notes and constraints

- `required_smarts` and `must_not_have_smarts` are present in the proto but are
  not implemented. For structural queries, it is usually better to combine this
  tool with `tsubstructure`.
- `-c` cannot be combined with `max_chiral`; removing chirality and filtering on
  chirality are contradictory.
- Utility mode changes output format unless `-u` is used.
