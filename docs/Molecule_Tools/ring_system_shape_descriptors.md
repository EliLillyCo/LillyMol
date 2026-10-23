# ring_system_shape_descriptors

## Overview

`ring_system_shape_descriptors` identifies molecules with a rod-like arrangement
of ring systems and substituents. The tool writes a descriptor table for every
input molecule and can optionally write SMILES files containing the molecules
that pass or fail the rod-like selection rule.

The immediate use case is selecting molecules whose ring systems and linkers
support an elongated shape rather than a compact, highly branched, or angular
shape.

The core ring-system calculation asks whether the exit points from each ring
system are as far apart as possible within that ring system. For example, a
para-disubstituted benzene is rod-like, while ortho- and meta-disubstituted
benzenes have a rod deficit.

If only selected ring atoms should contribute exit vectors, label those atoms
with an isotope and use `-I <iso>`.

The molecular-level decision also excludes substantial branching outside ring
systems and branching at ring atoms.

## Quick start

Compute descriptors:

```shell
ring_system_shape_descriptors input.smi > shape.txt
```

Write rod-like molecules while also keeping the descriptor table:

```shell
ring_system_shape_descriptors -R rodlike.smi input.smi > shape.txt
```

Write both rod-like and non-rod-like molecules:

```shell
ring_system_shape_descriptors \
  -R rodlike.smi \
  -N non_rodlike.smi \
  input.smi > shape.txt
```

The `-R` and `-N` files contain normal two-column SMILES output:

```text
smiles identifier
```

## Descriptor columns

The output is a space-separated descriptor table by default. Use `-o tab` for
tab-separated output.

| Column | Meaning |
| --- | --- |
| `Name` | Molecule identifier. |
| `nringsys` | Number of ring systems. |
| `nringsys_terminal` | Ring systems with zero or one meaningful exit point. These do not violate rod-likeness. |
| `nringsys_applicable` | Ring systems with exactly two meaningful exit points. |
| `nringsys_rodlike` | Applicable ring systems whose exit points are maximally separated. |
| `nringsys_not_rodlike` | Applicable ring systems whose exit points are not maximally separated. |
| `nringsys_multisub` | Ring systems with more than two meaningful exit points. |
| `nringsys_invalid` | Ring systems that could not be analysed. |
| `non_ring_branch_count` | Count of substantial branch points outside ring systems. |
| `ring_atom_branch_count` | Count of ring atoms with two or more substantial outside connections. With `-I`, only labelled ring atoms are counted. |
| `rodlike_molecule` | Final molecule-level rod-like classification, `1` or `0`. |
| `rod_deficit_min` | Minimum deficit among applicable non-terminal ring systems. |
| `rod_deficit_max` | Maximum deficit among applicable non-terminal ring systems. |
| `rod_deficit_ave` | Average deficit among applicable non-terminal ring systems. |

## Rod-like molecule rule

A molecule is classified as rod-like when all of the following are true:

- the molecule contains at least one ring system;
- no applicable ring system is non-rod-like;
- no ring system is multi-substituted;
- no ring system is invalid;
- `non_ring_branch_count` is zero;
- `ring_atom_branch_count` is zero.

In descriptor-table terms, the most useful selection condition is generally:

```text
rodlike_molecule == 1 && non_ring_branch_count == 0
```

The second condition is currently redundant because `rodlike_molecule` already
requires `non_ring_branch_count == 0`, but it documents the intended operational
selection.

## Ring-system geometry

For each ring system, the tool identifies meaningful exit points: ring-system
atoms bonded to atoms outside that ring system. Terminal exocyclic multiply
bonded atoms, such as carbonyl oxygens, are ignored as exit points.

With `-I <iso>`, only ring atoms with isotope `<iso>` are eligible to define
exit points. Substituents attached to unlabelled ring atoms, or to ring atoms
with a different isotope, are ignored for ring-system geometry. The total ring
system count is unchanged; a ring system with no labelled exit-vector atoms is
therefore treated as terminal for the shape classification.

For a ring system with exactly two exit points, the tool asks whether each exit
point is among the farthest atoms from the other exit point when distances are
computed within that ring system. If so, the ring system is rod-like.

Examples:

| Pattern | Classification |
| --- | --- |
| Terminal phenyl | terminal, does not violate rod-likeness |
| Ortho-disubstituted benzene | not rod-like |
| Meta-disubstituted benzene | not rod-like |
| Para-disubstituted benzene | rod-like |
| Ring with more than two exit points | multi-substituted |

The rod deficit is the number of ring-system bonds by which the observed exit
separation falls short of the best available separation.

## Branching descriptors

The ring geometry calculation alone is not sufficient to identify globally
rod-like molecules. A ring may have optimal exit vectors while a substituent
branches immediately at a ring atom, or a linker outside the ring system may
branch substantially.

`non_ring_branch_count` counts substantial branch points outside ring systems.
Terminal single-atom decorations such as methyl, fluoro, hydroxy, or amino
groups are not considered substantial. Atoms with terminal multiply bonded
neighbours, such as amide carbonyl carbons or sulfone sulfurs, are also not
considered branch points.

`ring_atom_branch_count` counts ring atoms with two or more substantial outside
connections. This catches cases such as an aliphatic ring atom bearing two
larger substituents. The ring-system geometry may still be optimal, but the
molecule is not considered rod-like.

When `-I <iso>` is specified, `ring_atom_branch_count` is also restricted to
ring atoms with isotope `<iso>`. `non_ring_branch_count` is not isotope-filtered.

## Isotope-restricted exit vectors

Use `-I <iso>` when only a subset of ring atoms should be considered as possible
exit-vector atoms.

```shell
ring_system_shape_descriptors -I 1 input.smi > shape.txt
```

This is useful when upstream processing has labelled the ring or ring-system
atoms of interest. If every relevant ring atom has isotope `1`, the calculation
is restricted to those labelled ring systems. If only selected ring atoms have
isotope `1`, the isotope acts as a precise inclusion mask for exit vectors and
ring-atom branch points. Unlabelled ring systems are still counted in `nringsys`,
but they contribute no labelled exit vectors.

## Ignoring terminal ring decorations

The `-x` option removes terminal single-atom ring substituents before the
ring-system geometry calculation when the same ring system also has at least one
larger attachment.

```shell
ring_system_shape_descriptors -x input.smi > shape.txt
```

This is useful for ignoring small decorations such as `F`, `OH`, `NH2`, or
`CH3` when a larger substituent defines the relevant shape of the molecule.

The branching descriptors are computed before this optional removal, so `-x`
only affects the ring-system geometry part of the calculation.

## Writing selected molecules

`-R` writes molecules that satisfy the rod-like selection rule:

```shell
ring_system_shape_descriptors -R rodlike.smi input.smi > shape.txt
```

`-N` writes molecules that do not satisfy the rule:

```shell
ring_system_shape_descriptors -N non_rodlike.smi input.smi > shape.txt
```

Both may be specified:

```shell
ring_system_shape_descriptors \
  -R rodlike.smi \
  -N non_rodlike.smi \
  input.smi > shape.txt
```

The output molecules are written as the original input SMILES and identifier.
This matters because the calculation may make temporary structural changes, for
example when `-x` is active.

## Other options

Common LillyMol molecule input options are supported:

| Option | Meaning |
| --- | --- |
| `-i <type>` | Input type. Usually inferred from filename. |
| `-g ...` | Chemical standardisation. |
| `-l` | Reduce to largest fragment. |
| `-c` | Remove chirality. |
| `-A ...` | Aromaticity options. |
| `-E ...` | Element options. |
| `-I <iso>` | Only ring atoms with isotope `<iso>` can define exit vectors or ring-atom branch points. |
| `-o <sep>` | Output separator. Recognised names include `tab` and `space`. |
| `-v` | Verbose progress and summary reporting. |

With `-v`, the tool reports summary statistics including the number of molecules
classified as rod-like and the number with non-ring or ring-atom branch points.
