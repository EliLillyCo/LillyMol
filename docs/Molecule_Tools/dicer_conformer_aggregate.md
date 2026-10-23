# dicer_conformer_aggregate

`dicer_conformer_aggregate` reads per-molecule output from `dicer -I geom`,
groups equivalent fragments, aligns their 3D geometries, and removes duplicate
conformations.

The initial implementation is intended primarily for linkers. It retains
fragments with exactly two non-hydrogen attachment atoms by default.

## Example

```shell
dicer -B serialized_proto -S diced.tfdata -I 1 -I geom conformers.sdf
dicer_conformer_aggregate -i tfdata -o tfdata -S linkers.tfdata diced.tfdata
```

Textproto input and output are useful for inspection and testing:

```shell
dicer -B proto -I 1 -I geom conformers.sdf > diced.textproto
dicer_conformer_aggregate diced.textproto > linkers.textproto
```

Each output `dicer_conformers::FragmentConformerSet` contains the normalized
fragment unique SMILES, an augmented topology containing its external
attachment atoms, atom roles, occurrence counts, and repeated packed xyz
coordinate arrays. Accepted conformers are placed in the coordinate frame of
the first conformer.

## Duplicate detection

For each fragment topology, the tool performs one self-substructure search and
retains role-preserving symmetry embeddings. Incoming geometries are tried with
each embedding. Fragment heavy atoms define the rigid alignment and internal
RMSD. The same transformation is applied to attachment atoms, whose maximum
displacement is evaluated separately.

A geometry is a duplicate when both of these conditions hold:

- fragment heavy-atom RMSD is no greater than `-r` (default 0.25);
- maximum attachment-atom displacement is no greater than `-a` (default 0.50).

The tool requires at least three non-collinear fragment heavy atoms. Alignment
uses rotation and translation only; molecular geometries are never scaled.

## Explicit hydrogens

All explicit hydrogen atoms are removed. External attachment records whose
external atom is hydrogen are also ignored. Hydrogens therefore do not affect
the aggregation key, attachment count, symmetry embeddings, alignment, or
duplicate calculation.

## Options

- `-i textproto|tfdata`: input format; default `textproto`.
- `-o textproto|tfdata`: output format; default `textproto`.
- `-S file`: output file. Textproto defaults to standard output; TFDataRecord
  output requires this option.
- `-m n`, `-M n`: minimum and maximum attachment counts; both default to 2.
- `-r distance`: fragment heavy-atom RMSD tolerance.
- `-a distance`: maximum attachment-atom displacement.
- `-e n`: maximum number of self-symmetry embeddings; default 1000.
- `-v`: report processing statistics.

Textproto input is expected to contain one complete `DicedMolecule` message per
line, matching the output produced by `dicer -B proto`.
