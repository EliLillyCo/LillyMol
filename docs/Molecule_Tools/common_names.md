# common_names

## Contents

- [Overview](#overview)
- [Quick start](#quick-start)
- [How grouping works](#how-grouping-works)
- [Structural equivalence](#structural-equivalence)
- [Representative structures](#representative-structures)
- [Output modes](#output-modes)
- [Two-pass and single-pass operation](#two-pass-and-single-pass-operation)
- [Large datasets and parallel processing](#large-datasets-and-parallel-processing)
- [Option reference](#option-reference)
- [Typical workflows](#typical-workflows)
- [Algorithm and limitations](#algorithm-and-limitations)

## Overview

`common_names` groups structurally equivalent molecules and retains the names of
all molecules in each group. It is related to
[`unique_molecules`](unique_molecules.md), but the two tools answer different
questions:

- `unique_molecules` asks, “Have I seen this structure before?” It writes the
  first representative and discards later duplicates.
- `common_names` asks, “Which input records represent the same structure?” It
  writes one representative and combines the names of every member of the
  group.

Use `unique_molecules` when the objective is simply to remove duplicates. Use
`common_names` when duplicate provenance, identifiers, or group sizes are
important.

For example, given:

```text
C methane
C carbon
CC ethane
```

running:

```shell
common_names -S common input.smi
```

creates `common.smi`:

```text
C methane:carbon
CC ethane
```

The first structure represents the group, and its name contains the identifiers
of both equivalent input records.

## Quick start

Group molecules using their normal LillyMol unique SMILES:

```shell
common_names -S common input.smi
```

Ignore chirality and directional-bond information while grouping:

```shell
common_names -c -x -S common input.smi
```

Standardize molecules and compare only their largest fragments:

```shell
common_names -g all -l -S common input.smi
```

Write the first identifier and group size instead of every identifier:

```shell
common_names -y -S common input.smi
```

For normal two-pass operation, `-S` is required unless textproto output is
selected. Add `-v` for processing statistics.

## How grouping works

For every input molecule, `common_names`:

1. Applies the requested comparison transformations.
2. Generates a LillyMol unique SMILES comparison key.
3. Looks up that key in an in-memory hash.
4. Creates a new group or adds the molecule name to the existing group.
5. Writes one representative molecule for each key after grouping is complete.

The default separator between accumulated names is `:`. Change it with `-D`.
For example:

```shell
common_names -D '|' -S common input.smi
```

would produce a name such as:

```text
C methane|carbon
```

Names are retained in input order within each group.

## Structural equivalence

By default, equality is based on LillyMol unique SMILES. Options can deliberately
remove or transform structural information before the comparison key is formed.

| Option | Comparison behavior |
| --- | --- |
| none | normal unique SMILES |
| `-c` | ignore chirality |
| `-x` | ignore directional cis/trans bond information |
| `-l` | compare only the largest fragment |
| `-I` | ignore isotopic labels |
| `-T E1=E2` | transform element `E1` to `E2` |
| `-g ...` | apply chemical standardization |
| `-a` | compare molecular graph forms |
| `-a -a` | compare graph forms plus total implicit-hydrogen count |

These options define the groups; they do not imply that the original records
were identical in every chemical detail.

### Chirality and directional bonds

Use `-c` when enantiomers should belong to the same group. Use `-x` when
cis/trans or directional-bond information should be ignored. Directional-bond
annotations in source data are often incomplete or inconsistent, so `-x` is
frequently appropriate when stereochemical provenance is not required.

### Fragments, isotopes, and standardization

`-l` removes all but the largest fragment before generating the comparison key.
This is commonly used to group salt forms by their parent structure. `-I`
removes isotope distinctions.

`-g all` applies LillyMol chemical standardization. This is particularly useful
when records come from several sources that may represent functional groups or
charges differently.

### Element transformations

`-T` can make selected elements equivalent. For example:

```shell
common_names -T I=Cl -T Br=Cl -S common input.smi
```

groups iodine, bromine, and chlorine forms according to the transformed
comparison structure. Fluorine remains distinct.

### Graph comparison

`-a` converts molecules to graph form before comparison: bond orders, charges,
chirality, and isotopes are removed. This is a broad equivalence relation; for
example, benzene and cyclohexane have the same carbon-ring graph.

Specify `-a` twice to append the starting molecule’s total implicit-hydrogen
count to the graph key:

```shell
common_names -a -a -S common input.smi
```

This distinguishes benzene from cyclohexane while retaining graph-based
comparison.

## Representative structures

Each group is represented by the first matching molecule encountered in the
input. The structure written is the original representation of that first
molecule, even when a transformed structure was used as the comparison key.

For example, given:

```text
BrC bromo
ClC chloro
FC fluoro
IC iodo
```

and:

```shell
common_names -T I=Cl -T Br=Cl -S common input.smi
```

the output includes:

```text
BrC bromo:chloro:iodo
FC fluoro
```

The comparison key for the first group uses transformed chlorine, but the output
SMILES remains `BrC` because that was the first representative. Input order can
therefore affect the displayed structure and the order of names, but not group
membership.

## Output modes

### Concatenated names

The default output writes one molecule per group and concatenates all names with
the separator specified by `-D`, or `:` by default:

```shell
common_names -D ':' -S common input.smi
```

`-S common` writes `common.smi` unless another output type is selected with
`-o`.

In normal two-pass mode, do not choose an output stem that resolves to the same
file as an input file. The output is opened before the second input pass, so
using the input basename as the output stem can overwrite the file that still
needs to be reread.

### First name and count

`-y` retains only the first name and the number of instances:

```text
C methane 2
```

This can save substantial memory when groups contain many names. It is available
in normal two-pass mode, not with `-f`.

### Concatenated names and count

`-X num` appends a separate group-count column while retaining the concatenated
names:

```shell
common_names -X num -S common input.smi
```

### Textproto output

`-X proto=<file>` writes one `common_names::CommonNames` textproto per group:

```shell
common_names -X proto=groups.textproto input.smi
```

A record resembles:

```textproto
smiles: "[N+]1(=CN(C)C=C1)CCCC"
key: "CCCC[n+]1c[n](C)cc1"
id: "CHEMBL3182180"
id: "CHEMBL3184676"
id: "CHEMBL3559949"
```

`smiles` is the first representative, `key` is the transformed comparison key,
and repeated `id` fields hold the names in that group. Textproto output cannot
be combined with regular `-S` output.

## Two-pass and single-pass operation

### Normal two-pass mode

By default, `common_names` reads its input twice:

1. The first pass generates comparison keys and accumulates names or counts.
2. The second pass writes the first representative of each group in input order.

To avoid recomputing and retaining full molecules, this mode stores one
comparison key per input record. It therefore needs to know the number of input
molecules before allocating its arrays. Without `-s`, the tool counts the input
records before grouping them. With `-s <size>`, the supplied value is used as the
allocation size and must not be smaller than the actual number of molecules.

Two-pass mode requires seekable input files and cannot consume a non-replayable
stdin stream.

### Single-pass mode

`-f` reads the input once and retains enough information to write SMILES output
without rereading the source:

```shell
generate_molecules | common_names -f -g all -l -c -S - -
```

Single-pass mode:

- supports stdin;
- writes SMILES only;
- ignores `-o` and `-s`;
- uses more memory because it retains the first input SMILES and comparison key
  information needed for output.

Use it for pipelines or when rereading the source is impossible. For large,
seekable files, normal two-pass mode generally has the lower memory cost.

## Large datasets and parallel processing

`common_names` retains one hash entry per distinct comparison key and either all
names or a count for each group. Memory therefore grows with both structural
diversity and, unless `-y` is used, the total identifier volume.

For large single-file datasets, `unique_molecules_parallel.sh` can run
`common_names` on independent formula partitions:

```shell
unique_molecules_parallel.sh \
  -thr 16 -exe common_names.sh -S common \
  -c -x input.smi
```

The wrapper first uses `msort_parallel` to sort and split the input by molecular
formula, runs one `common_names` process per partition, and concatenates the
partition outputs. Despite its name, the wrapper can run either
`unique_molecules` or `common_names` through `-exe`.

Important constraints are:

- only one input file is accepted;
- stdin is not supported;
- output order is not guaranteed to match input order;
- temporary storage can approach several copies of the input size;
- `-nj` leaves partition outputs separate rather than concatenating them.

Most importantly, formula partitioning is correct only when molecules considered
equivalent by `common_names` necessarily had the same **original** molecular
formula. Options such as `-l`, `-T`, `-a`, or formula-changing standardizations can make
molecules with different original formulae equivalent; such records may be sent
to different partitions and will not be grouped together. Use serial
`common_names`, or preprocess the input before partitioning, when comparison
transformations can cross original-formula boundaries.

For additional wrapper details, see
[`unique_molecules_parallel`](unique_molecules_parallel.md).

## Option reference

### Comparison controls

- `-c`: ignore chirality.
- `-x`: ignore directional cis/trans bond information.
- `-l`: retain only the largest fragment for comparison.
- `-I`: remove isotopes for comparison.
- `-T E1=E2`: transform elements before comparison; may be repeated.
- `-g ...`: apply chemical standardization.
- `-a`: compare graph forms; specify twice to include implicit-hydrogen count.

### Group and output controls

- `-D <separator>`: separator used between accumulated names; default `:`.
- `-S <stem>`: regular output stem; required in normal mode.
- `-o <type>`: output molecule format in normal mode.
- `-y`: write first name and instance count instead of all names.
- `-X num`: append the group count while retaining concatenated names.
- `-X proto=<file>`: write grouped records as textproto instead of molecule
  output.
- `-f`: use single-pass, SMILES-only operation.

### Input sizing and progress

- `-s <size>`: preallocate for at most `<size>` input molecules in normal mode.
  The value must not be smaller than the actual input size.
- `-r <count>`: report progress every `<count>` molecules.
- `-p <file>`: seed grouping with a previously collected molecule file.
- `-i <type>`: specify the input molecule type.
- `-v`: increase diagnostic output.

Standard LillyMol `-A`, `-E`, and `-K` controls configure aromaticity, elements,
and SMILES generation.

## Typical workflows

### Find exact structural groups

```shell
common_names -S common input.smi
```

### Group standardized parent structures

```shell
common_names -g all -l -S common input.smi
```

### Ignore stereochemistry

```shell
common_names -c -x -S common input.smi
```

### Compare graph forms while retaining hydrogen count

```shell
common_names -a -a -S common input.smi
```

### Write names and group counts

```shell
common_names -X num -S common input.smi
```

### Reduce identifier memory

```shell
common_names -y -S common input.smi
```

### Generate a machine-readable group report

```shell
common_names -X proto=groups.textproto input.smi
```

## Algorithm and limitations

Grouping uses an in-memory hash from comparison key to accumulated names or
count. Normal mode also stores the comparison key associated with every input
record so it can identify first representatives during the second pass.

Consequences include:

- memory use scales with the number of records, distinct structures, and retained
  identifier text;
- normal mode rereads every input file and therefore requires unchanged,
  seekable inputs;
- the first input molecule determines the output representation of a group;
- names containing the selected separator can be ambiguous in concatenated
  output; textproto output is preferable when identifiers may contain delimiter
  characters;
- graph comparison and aggressive transformations create deliberately broad
  equivalence classes and should be chosen with chemical intent.

Use `unique_molecules` instead when duplicate names and group provenance are not
needed.
