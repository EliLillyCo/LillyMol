# chemotypes

## Overview

`chemotypes` extracts query-defined ring-system chemotypes from molecules. A
substructure query identifies the seed ring system, and the output molecule is
reduced to that ring system plus any optional neighbouring ring systems,
linkers, exocyclic double-bonded atoms, or attachment-point annotations.

The intended use is grouping molecules by a chemically meaningful ring-defined
core. Chemists often describe a series by a chemotype such as a piperidine,
benzimidazole, or substituted fused ring system. This tool makes that
definition explicit by combining a query with deterministic rules for expanding
from the matched seed ring system.

## Quick Start

Extract the ring system containing the first atom matched by a SMARTS:

```shell
chemotypes -s '[N]' -z i -i smi input.smi > chemotypes.smi
```

For example, from:

```text
CN1CCCCC1 mol1
CC miss
```

the command writes:

```text
N1CCCCC1 mol1
```

By default, a molecule that does not match any query is an error. Use `-z i`
or `-z ignore` to ignore nonmatching molecules.

## Query Matching

The seed must be specified with one or more queries:

```shell
chemotypes -q seed.qry input.smi > chemotypes.smi
chemotypes -s '[nH]1cccc1' input.smi > chemotypes.smi
```

Queries are examined in command-line order. The first matching query is used.
By default, that query must identify exactly one ring system across all of its
embeddings. This means a broad query such as `c` is acceptable for a molecule
with one aromatic ring system, but ambiguous for biphenyl because the embeddings
touch two different ring systems.

Within the accepted query, the first embedding that touches the identified ring
system is used. The first matched atom in that embedding that belongs to the
ring system defines the seed atom.

Use `-z f` or `-z first` to keep the historical behaviour of using the first
embedding of the first matching query, even if the query touches multiple ring
systems.

The query should normally include ring constraints where appropriate. If the
first matching query has no matched atom in a ring system, processing stops for
that molecule.

## Adjacent Ring Systems

Use `-n <n>` to include directly adjacent ring systems.

```shell
chemotypes -s '[c]' -n 1 -i smi input.smi > chemotypes.smi
```

Ring systems are considered directly adjacent when they can be reached from the
seed ring system through a path whose internal atoms are not ring atoms. Ring
systems beyond another intervening ring system are not considered directly
adjacent.

The adjacent systems are ranked by distance from the seed ring system, with
query atom order used to make deterministic choices in ambiguous cases. If a
molecule has more adjacent systems than requested, only the best ranked systems
are included.

Use `-Y larger` to prefer the same-distance candidate whose candidate side
contains more atoms. The size includes linker atoms from the seed ring system to
the candidate ring system and the atoms reachable through the bond entering that
candidate ring system. Distance remains the primary ranking criterion;
`-Y larger` only resolves same-distance choices.

Use `-t` with `-n` to include all systems tied at the cutoff distance:

```shell
chemotypes -s '[c]' -n 1 -t -i smi input.smi > chemotypes.smi
```

In a case where two phenyl rings are directly attached to the seed ring at the
same distance, `-n 1` includes one of them, while `-n 1 -t` includes both.

Use `-D <dist>` to limit adjacent ring systems by maximum bond distance from
the seed ring system.

## Linkers and Exocyclic Atoms

When adjacent ring systems are included, linker atoms between selected ring
systems are included by default. The linker is the non-ring connected component
that connects the selected ring systems, not just one shortest path.

Terminal exocyclic atoms double bonded to retained atoms are also included by
default. This preserves common ring-system features such as carbonyl oxygen
atoms.

## Attachment Points

The default output contains only the retained chemotype core and does not encode
how that core was attached to the rest of the parent molecule.

Use `-I <iso>` to label retained atoms at bonds crossing the final chemotype
boundary:

```shell
chemotypes -s '[N]' -I 99 -i smi input.smi > labelled.smi
```

This labels the retained endpoint of every removed bond with the fixed isotope
value, including cuts beyond atoms retained by `-u`. `-I` cannot be used with
atom typing.

Use `-C <stem>` to also write the complement to a separate file:

```shell
chemotypes -s '[N]' -I 99 -C complement -i smi input.smi > chemotypes.smi
```

The complement contains exactly the excluded atoms. Both endpoints of each cut
bond receive isotope 99, so the attachment atoms are marked in both outputs.
The complement uses the `-o` output types (SMILES by default), preserves the
parent name, and keeps disconnected pieces together in one molecule record.
No record is written for an empty complement. `-C` requires a separate file
stem and cannot be used with fingerprint output. A fixed isotope marks sites;
it does not uniquely pair multiple connections.

Use `-u` to retain one-hop atoms attached to retained ring atoms without
applying atom-type labels:

```shell
chemotypes -s '[N]' -u -i smi input.smi > attached.smi
```

When `-P` is specified, atom typing is applied before the molecule is reduced.
The tool automatically retains one-hop atoms attached to the chemotype core and
keeps isotopic atom-type labels on non-terminal attachment atoms. Isotopes on
the core chemotype atoms and on atoms that were singly connected in the parent
molecule are cleared. This can distinguish, for example, the same ring system
attached through carbon or nitrogen.

```shell
chemotypes -s '[N]' -P UST:ARY -i smi input.smi > typed.smi
```

Use `-x` with `-u` to ignore singly connected attached atoms when retaining
untyped attachment atoms:

```shell
chemotypes -s '[N]' -u -x -i smi input.smi > attached.smi
```

## Fingerprint Output

Use `-J` to write a fingerprint for the chemotype atom set instead of writing
chemotype molecules. Fingerprint mode writes TDT records to stdout and disables
the normal molecule output and summary options. If no `-P` atom type is
specified, `UST:ARY` is used.

A tag beginning with `FP` generates a fixed-width linear fingerprint using the
standard LillyMol linear fingerprint generator:

```shell
chemotypes -s '[c]' -J FPCTYPE -i smi input.smi > chemotype.gfp
```

A tag beginning with `NC` generates a non-colliding counted EC fingerprint:

```shell
chemotypes -s '[c]' -J NCCTYPE -i smi input.smi > chemotype.gfp
```

For `NC` fingerprints, a trailing digit on the tag sets the EC radius. If no
trailing digit is present, radius 3 is used.

```shell
chemotypes -s '[c]' -J NCCTYPE2 -i smi input.smi > chemotype.gfp
```

By default, only atoms in the chemotype are fingerprinted. Use
`-J EXPAND=<n>` to include atoms within `n` bonds of that selected atom set:

```shell
chemotypes -s '[c]' -J NCCTYPE2 -J EXPAND=1 -i smi input.smi > chemotype.gfp
```

Use `-J INVERT`, or the synonym `-J OUTSIDE`, to fingerprint the atoms outside
the chemotype instead:

```shell
chemotypes -s '[c]' -J FPCTYPE -J OUTSIDE -i smi input.smi > outside.gfp
```

With `-z i` or `-z ignore`, a molecule that does not match the query writes an
empty fingerprint in a normal TDT record. Without `-z i`, nonmatching molecules
remain fatal.

Use `-f` to insert the chemotype fingerprint into an existing TDT stream. The
input is copied to stdout, and the new fingerprint is written immediately after
each `$SMI<...>` record.

```shell
chemotypes -s '[c]' -J FPCTYPE -f input.gfp > with_chemotype.gfp
```

## Parent Output

Use `-p <text>` to write the parent molecule before each generated chemotype:

```shell
chemotypes -s '[N]' -p PARENT -z i -i smi input.smi > parent_and_chemotype.smi
```

The value is appended to the parent molecule name as a separate token, which
makes parent records easier to distinguish from chemotype records when reviewing
SMILES output. Use `-p .`, `-p def`, or `-p default` to write the parent without
adding an annotation.

This is mainly useful for debugging, review, and explaining why a chemotype was
generated. The parent is written only for molecules that successfully generate a
chemotype.


## Chemotype Hashes

Use `-H def` to append a deterministic fixed-width alphanumeric hash to each
written chemotype record:

```shell
chemotypes -s '[N]' -z i -i smi -H def input.smi > chemotypes.smi
```

The hash is computed from the unique SMILES of the reduced chemotype, so the
same chemotype receives the same identifier regardless of which parent molecule
or input file generated it. The default hash is a 26-character uppercase base32
encoding derived from an MD5 digest. Use `-H width=<n>` to request a shorter
fixed width when a compact identifier is more important than collision margin.

The tool checks for collisions among chemotypes seen in the current run. A
collision, where two different unique chemotype SMILES map to the same hash, is
a fatal error.

Use `-H xref=<fname>` to also write a cross-reference file:

```shell
chemotypes -s '[N]' -z i -i smi -H xref=chemotype_hashes.txt input.smi > chemotypes.smi
```

The cross-reference file contains the hash, unique chemotype SMILES, count,
and first parent molecule name. The parent name is written last so the first
three fields remain easy to parse when parent names contain spaces.

## Summary Output

Use `-F <fname>` to accumulate generated chemotypes as
`dicer_data::DicerFragment` textprotos:

```shell
chemotypes -s '[N]' -z i -i smi -F summary.textproto input.smi > chemotypes.smi
```

The summary is keyed by the unique SMILES of each generated chemotype. The first
occurrence creates a proto with:

| Field | Meaning |
| --- | --- |
| `smi` | Unique SMILES of the chemotype. |
| `par` | Name of the first parent molecule that generated the chemotype. |
| `nat` | Number of atoms in the chemotype. |
| `n` | Number of parent molecules that generated this chemotype. |
| `iso` | Present when isotopic labels encode attachment information. |

For example, if two input molecules reduce to the same piperidine chemotype,
the summary contains:

```text
smi: "N1CCCCC1" par: "mol1" nat: 6 n: 2
```

When `-P` is active, `iso` is written as `ATYPE`. When `-I` is active, `iso` is
written as `ATT`.

## Filtering by Ring Count

Use `-r <n>` to process only molecules with at least `n` rings:

```shell
chemotypes -s '[c]' -n 1 -r 2 -z i -i smi input.smi > chemotypes.smi
```

This is useful when the requested chemotype expansion requires nearby rings and
single-ring molecules should simply be skipped.

## Options

| Option | Meaning |
| --- | --- |
| `-q <query>` | LillyMol query file defining the chemotype seed. May be repeated. |
| `-s <smarts>` | SMARTS defining the chemotype seed. May be repeated. |
| `-n <n>` | Include `n` directly adjacent ring systems. Default `0`. |
| `-t` | With `-n`, include all adjacent ring systems tied at the cutoff distance. |
| `-Y larger` | When selecting adjacent ring systems with `-n`, prefer the same-distance candidate whose candidate side contains more atoms. Synonyms: `largest`, `size`. |
| `-D <dist>` | Maximum bond distance to an adjacent ring system. |
| `-r <n>` | Minimum number of rings required for processing. |
| `-u` | Include one-hop atoms attached to retained ring atoms. |
| `-x` | With `-u`, ignore singly connected attached atoms. |
| `-C <stem>` | Write excluded atoms to a separate output stem using the `-o` types. |
| `-I <iso>` | Label retained boundary atoms, and complement endpoints, with isotope `iso`. Incompatible with `-P`. |
| `-P <atype>` | Atom typing specification for retained non-terminal attachment atoms; also controls fingerprint atom types when `-J` is active. |
| `-J FP<tag>` | Write a fixed-width linear fingerprint of the chemotype atom set. |
| `-J NC<tag>` | Write a non-colliding counted EC fingerprint of the chemotype atom set; trailing digits set the EC radius. |
| `-J EXPAND=<n>` | Expand the fingerprinted atom set by `n` bonds. |
| `-J INVERT`, `-J OUTSIDE` | Fingerprint atoms outside the chemotype. |
| `-f` | Work as a TDT filter, inserting the `-J` fingerprint after each `$SMI<...>` record. |
| `-p <text>` | Write the parent molecule before each generated chemotype; append `<text>` unless it is `.`, `def`, or `default`. |
| `-F <fname>` | Write accumulated `dicer_data::DicerFragment` textproto summary. |
| `-H def` | Append a deterministic fixed-width alphanumeric hash of the unique chemotype SMILES. |
| `-H xref=<fname>` | Also write a hash cross-reference file. |
| `-H width=<n>` | Set hash width. Default and maximum are 26 characters. |
| `-z i`, `-z ignore` | Ignore molecules that do not match any query. |
| `-z f`, `-z first` | Use the first embedding when a query matches multiple ring systems. |
| `-S <stem>` | Output file stem. Default stdout. |
| `-i <type>` | Input type. Usually inferred from filename. |
| `-o <type>` | Output type. Default SMILES. |
| `-g ...` | Chemical standardisation. |
| `-l` | Reduce to largest fragment. |
| `-c` | Remove chirality. |
| `-A ...` | Aromaticity options. |
| `-E ...` | Element options. |
| `-v` | Verbose progress and summary reporting. |

