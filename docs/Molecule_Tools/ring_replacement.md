# Ring Replacement

## Contents

- [Overview](#overview)
- [Quick start](#quick-start)
- [How ring replacement works](#how-ring-replacement-works)
- [Worked example and output](#worked-example-and-output)
- [Choosing breadth versus confidence](#choosing-breadth-versus-confidence)
- [Using `ring_replacement`](#using-ring_replacement)
- [Atom typing attachment context](#atom-typing-attachment-context)
- [Controlling molecular-formula difference](#controlling-molecular-formula-difference)
- [Preserving 3D coordinates](#preserving-3d-coordinates)
- [Python](#python)
- [Using `ring_replacement_inexact`](#using-ring_replacement_inexact)
- [Building replacement-ring collections](#building-replacement-ring-collections)
  - [Extracting rings](#extracting-rings-with-ring_extraction)
  - [Replacement-file naming](#replacement-file-naming)
  - [Administrator wrapper](#administrator-wrapper)
  - [File format](#file-format)
  - [Spiro systems and exocyclic double bonds](#spiro-systems-and-exocyclic-double-bonds)
  - [Extracting several source collections](#extracting-several-source-collections)
  - [Aggregating collections](#aggregating-collections)

## Overview

`ring_replacement` generates analogues by replacing a ring or ring system in an
input molecule with rings observed in existing molecular collections. This is a
form of **scaffold hopping**: the central ring framework changes while the
surrounding substituents are retained.

For example, an indole ring system might be replaced by a benzofuran,
benzothiophene, indazole, or another fused five- and six-membered aromatic
heterocycle. Bonds from the original ring system to the rest of the molecule are
retained, preserving the substituents attached to the original scaffold.

Unlike a similarity search, this operation proposes new structures rather than
retrieving existing molecules. The replacement collections record both the ring
structure and an exemplar count, allowing frequently observed rings to be
preferred over unusual or potentially erroneous structures.

The standard replacement operation preserves ring size, aromaticity, and
substitution pattern. A para-disubstituted six-membered aromatic ring is therefore
replaced by other para-disubstituted six-membered aromatic rings.

## Quick start

This section assumes that a collection of replacement rings is already
available. LillyMol includes replacement rings derived from ChEMBL under
`${LILLYMOL_HOME}/data/ring_replacement`. The following command replaces
six-membered aromatic rings using replacements with at least five precedents:

```shell
ring_replacement \
  -R "${LILLYMOL_HOME}/data/ring_replacement/rings_6a.smi" \
  -n 5 -u -I . input.smi > variants.smi
```

The important options are:

- `-R` specifies a replacement-ring collection produced by `ring_extraction`.
  It may be specified more than once.
- `-n 5` considers only rings observed in at least five source molecules.
- `-u` suppresses duplicate products.
- `-I .` removes the temporary isotopic labels used while assembling products.

Without `-s` or `-q`, every compatible ring in each input molecule can be
replaced. Use `-s` or `-q` when only particular rings should be considered:

```shell
ring_replacement \
  -R "${LILLYMOL_HOME}/data/ring_replacement/rings_6a.smi" \
  -s '[/IWfss1cr6]' -z i -n 5 -u -I . input.smi > variants.smi
```

Here the SMARTS restricts replacement to isolated six-membered aromatic ring
systems. `-z i` quietly skips input molecules that do not match the selection
query.

## How ring replacement works

A replacement-ring record contains two aligned representations: a query that
identifies a compatible ring in the input molecule and a replacement structure
whose atoms can be substituted one-for-one. This alignment allows the tool to
retain bonds from the ring to the rest of the parent molecule.

By default, ring size, aromaticity pattern, and substitution pattern are
preserved. For example, replacing a para-disubstituted six-membered aromatic
ring produces other para-disubstituted six-membered aromatic rings. The
replacement changes the ring framework while retaining the surrounding
substituents.

A typical run therefore consists of:

1. Choose replacement files for the ring or ring-system classes of interest.
2. Optionally use `-s` or `-q` to select which parent rings may be replaced.
3. Restrict candidate rings by precedent, formula difference, or structural
   filters.
4. Generate products and optionally filter them with `-Y` and `-N` queries.

The support value stored with each replacement is the number of source
molecules containing that ring. It is evidence that the ring has occurred in a
real molecular collection, not proof that every generated product is
synthetically accessible.

## Worked example and output

Given the input molecule:

```text
Clc1ccc(F)cc1 p-benzene
```

a basic replacement using the six-membered aromatic collection is:

```shell
ring_replacement \
  -R "${LILLYMOL_HOME}/data/ring_replacement/rings_6a.smi" \
  -u -w p-benzene.smi > variants.smi
```

Example output includes:

```text
Cl[1c]1cc[1c](F)cc1 p-benzene %% CHEMBL503634.6a 1496512
Cl[1c]1c[n][1c](F)cc1 p-benzene %% CHEMBL156037.6a 83606
Cl[1c]1[n]c[1c](F)c[n]1 p-benzene %% CHEMBL1171471.6a 19926
Cl[1c]1[n][n][1c](F)cc1 p-benzene %% CHEMBL600052.6a 10360
```

Each output record contains:

1. The product SMILES.
2. The name of the starting molecule.
3. `%%`, separating parent and replacement metadata.
4. The identifier of an exemplar containing the replacement ring. The `6a`
   suffix denotes a six-membered aromatic ring.
5. The number of source molecules containing that ring.

The isotopes in this output identify attachment atoms used during replacement.
Specify `-I .` to remove them. Symmetry can allow one replacement ring to produce
more than one distinct product, and `-u` removes products that are identical.

The support count is a useful confidence signal. Rings with many precedents are
more likely to represent established chemistry; rings with very low counts may
be valid rare structures, drawing errors, or speculative chemistry. Use `-n` to
set an appropriate minimum for the application.

## Choosing breadth versus confidence

Ring replacement offers a progression from conservative generation to broad
scaffold exploration. In this discussion, *lower risk* or *safer* means that a
product is more strongly supported by precedents and therefore more likely to be
chemically realistic and makable. It does not refer to laboratory or toxicological
safety, and no setting guarantees synthetic accessibility.

From most conservative to most exploratory:

1. **Atom-typed exact replacement.** The ring-system topology and observed
   substitution pattern must match, and the atom types of the attached
   substituents must also be compatible. A 1,3-disubstituted ring is replaced
   only by a 1,3-disubstituted ring observed with the same kinds of substituents.
   This generally produces the fewest, best-supported products.
2. **Default exact replacement.** The topology and observed substitution pattern
   are preserved, but the kinds of substituents may differ from those in the
   precedent. This produces more products with less contextual evidence.
3. **Exact replacement with `-d`.** The basic ring-system topology is preserved,
   but the replacement need not have been observed with the product's
   substitution pattern. This increases product count and the risk of proposing
   unrealistic substitution patterns.
4. **Inexact replacement with attachment constraints.**
   `ring_replacement_inexact` permits a different ring-system shape or size. Its
   `-e` option preserves available same-ring attachment relationships from the
   starting scaffold while the possible placements are enumerated.
5. **Unconstrained inexact replacement.** A different ring-system topology is
   inserted and the possible substituent placements are enumerated without
   preserving those initial attachment relationships. This gives the broadest
   scaffold hopping and the largest, highest-risk product set.

The appropriate level depends on the purpose. Conservative modes are suitable
when synthetic plausibility is paramount. Less constrained modes are useful when
exploration and structural novelty justify reviewing a larger number of weaker
proposals.

## Using `ring_replacement`

`ring_replacement` performs exact replacement: the basic topology of the parent
ring or ring system and the replacement ring is the same. Its default behavior
also requires the observed substitution pattern to be preserved.

### Replacement-ring input

- `-R <file>` reads a replacement-ring collection generated by
  `ring_extraction`. The option may be repeated to use multiple collections.
- `-n <count>` considers only replacement rings with at least this many source
  examples.
- `-P <atype>` applies atom typing to attachment points. The same atom typing
  must have been used when the replacement collection was generated.

### Selecting rings in the starting molecule

- `-s <SMARTS>` selects rings containing atoms matched by the SMARTS. It may be
  specified more than once.
- `-q <query>` provides queries using the standard LillyMol query syntax.
- `-z i` ignores input molecules that do not match any `-s` or `-q` query.
- `-z w` writes nonmatching input molecules unchanged. It is normally combined
  with `-z i`.
- `-d` stops requiring preservation of an observed substitution pattern. The
  ring-system topology is still preserved.

If neither `-s` nor `-q` is specified, every ring compatible with the supplied
replacement collections can be considered.

### Filtering replacement candidates

- `-F <query>` requires candidate replacement rings to match a query.
- `-D <query>` rejects candidate replacement rings that match a query.
- `-f <difference>` limits the molecular-formula difference between the parent
  ring and replacement ring.

`-F` and `-D` operate on the replacement rings before products are generated.
They are useful for requirements such as retaining or excluding particular ring
heteroatoms.

### Filtering generated products

- `-Y <query>` requires generated products to match a query.
- `-N <query>` rejects generated products that match a query.

These options operate on complete product molecules. The same filtering could be
performed by a subsequent `tsubstructure` invocation, but applying it during
replacement avoids writing unwanted products.

### Output control

- `-u` writes unique products only.
- `-p` writes each parent molecule as well as its generated products.
- `-w` orders products by replacement-ring precedent count.
- `-I .` removes isotopic attachment labels from products.
- `-I <n>` changes existing isotopes to `<n>`.
- `-B <file>` writes molecules for which no replacement was generated to a
  separate SMILES file.
- `-3` transfers coordinates from corresponding parent-ring atoms to replacement
  atoms and writes coordinates in LillyMol SMILES form.

### Input preprocessing and diagnostics

- `-c` removes chirality before replacement.
- `-l` reduces input molecules to their largest fragment.
- `-g ...` applies chemical standardization.
- `-A ...`, `-E ...`, and `-i ...` provide the standard LillyMol aromaticity,
  element, and input-type controls.
- `-v` increases diagnostic output and may be repeated.
- `-X help` lists miscellaneous controls.


## Atom typing attachment context

The most conservative exact replacement uses atom typing to preserve not only
the substitution pattern, but also the kinds of atoms attached to the ring. The
atom type stored at a ring attachment point is the type of the **exocyclic atom**
that was attached there in the source molecule; it is not the type of the ring
atom itself.

For example, without atom typing, a replacement ring known with a 1,3
substitution pattern can be used for any compatible 1,3-disubstituted parent.
With atom typing, a parent bearing oxygen and carbon substituents at those
positions matches only precedents observed with compatible oxygen and carbon
attachments. This produces fewer products, but gives stronger evidence for the
local chemical context of each proposed attachment.

Atom typing must be enabled with the identical `-P` specification during both
ring extraction and replacement:

```shell
ring_extraction -A 2 -S RINGS -P UST:CY -k -c -l -v chembl.smi

ring_replacement \
  -R RINGS_5a6a.smi -P UST:CY -n 5 -u -I . \
  input.smi > variants.smi
```

If the specifications differ, the isotopic atom-type labels in the replacement
collection will not correspond to the types calculated for the input molecule,
and expected replacements will not match.

`UST:CY` is a useful starting point:

- `C` records the number of connections to the exocyclic atom.
- `Y` records a compressed atomic number in which related heavy halogens are
  treated equivalently.

A more restrictive type such as `UST:ABCHY` additionally records aromaticity,
unsaturation, and hydrogen count. Greater specificity increases contextual
confidence but can sharply reduce coverage, particularly when the source
collection is small. Atom typing is described in
[atom typing](/docs/Molecule_Lib/atom_typing.md).

## Controlling molecular-formula difference

The `-f <difference>` option limits how far the replacement ring may diverge
from the ring being removed. This can be important for common rings and ring
systems, where a large replacement collection can otherwise generate very many
products.

A small threshold emphasizes **exploitation** around the starting molecule. This
is useful when the parent is already interesting—for example, a hit from a
biological screen—and the goal is to explore nearby scaffold variants without
changing composition substantially. A larger threshold permits more
**exploration** of chemically different scaffold forms, increasing both product
count and structural diversity.

The difference is calculated for the ring or ring system, not for the complete
molecule. LillyMol counts changes across atom categories, including aromatic and
aliphatic forms, hydrogen environments, and ring atoms, and sums the absolute
count differences. For example, changing benzene to pyridine has a difference of
three: one fewer aromatic carbon, one additional aromatic nitrogen, and one fewer
aromatic carbon bearing hydrogen.

The tool rejects a threshold of one; use a value of at least two. A typical
invocation that keeps pyridine replacements relatively close is:

```shell
ring_replacement \
  -R "${LILLYMOL_HOME}/data/ring_replacement/rings_6a.smi" \
  -s '[/IWfss1c]1ncccc1' -z i -f 4 -p -u -I . \
  input.smi > variants.smi
```

Here `-f 4` accepts replacement rings whose calculated difference is no greater
than four. It should be interpreted as a distance threshold, not simply as the
number of atoms changed. Run with several thresholds to determine an appropriate
balance between product count and scaffold diversity for a particular collection.

## Preserving 3D coordinates

If input molecules contain 3D coordinates, `-3` transfers the coordinates of the
corresponding parent-ring atoms to the replacement atoms. This can be useful when
exploring analogues of a docked or otherwise positioned molecule. It is a
coordinate transfer, not a geometry optimization; new products may still require
relaxation or minimization.

With `-3`, coordinates are written using LillyMol's coordinate-bearing SMILES
syntax. A record begins like this:

```text
O{{0.0021,-0.0041,0.002}}=[N+]{{-0.0158,1.3128,0.0093}}([O-]{{1.1157,1.9867,0.0021}})C{{-1.3076,2.0354,0.0193}}...
```

The ellipsis deliberately truncates the example, so the text shown above is not
a valid SMILES record. A complete coordinate-bearing output file can be converted
to SDF form with:

```shell
fileconv -i smi -o sdf variants.smi
```

## Python

The `lillymol_tools` module exposes the common exact-replacement controls through
the `RingReplacement` class. `LILLYMOL_HOME` must be set so that a replacement
filename such as `rings_6a.smi` can be resolved under
`${LILLYMOL_HOME}/data/ring_replacement`.

```python
from lillymol import MolFromSmiles
from lillymol_tools import RingReplacement

molecule = MolFromSmiles("Oc1ccc(OC)cc1 start")

replacement = RingReplacement()
replacement.set_ring_atom_smarts("c")
replacement.set_min_support_requirement(100)
replacement.set_unique_molecules_only(True)
replacement.set_remove_isotopes(True)
# Optional: broaden replacement to substitution patterns not observed in the
# replacement collection. This is the Python equivalent of command-line -d.
replacement.set_preserve_connectivity(False)

number_read = replacement.read_replacement_rings("rings_6a.smi")
if number_read == 0:
  raise RuntimeError("Cannot read replacement rings")

products = replacement.process(molecule)
for product in products:
  print(product.smiles(), product.name())
```

`set_ring_atom_smarts("c")` means that rings containing a matched aromatic
carbon may be replaced. Use a more selective SMARTS when only a particular ring
or ring environment should be processed.

Configuration that filters replacement records, especially
`set_min_support_requirement`, should be applied before
`read_replacement_rings`. The Python API currently reads one replacement-ring
collection per `RingReplacement` object. `number_replacement_rings()` reports how
many records were retained from that collection.

By default, `set_preserve_connectivity(True)` requires replacement records to
preserve an observed substitution pattern. Calling
`set_preserve_connectivity(False)` is equivalent to command-line `-d`: it admits
records marked `conn: false`, broadening the generated structures and potentially
producing less plausible substitution patterns. This setting determines which
records are retained while the replacement collection is read, so it must be set
before `read_replacement_rings()`. Calling it after the collection has been read
raises `RuntimeError`.

`process` returns a Python list of generated `Molecule` objects. The parent
molecule is **not** included, and an empty list means that no product was
generated. A generated product can still be structurally identical to the parent
if the collection contains the same ring framework.

When `set_unique_molecules_only(True)` is active, duplicate tracking persists
across calls to `process` on the same object. This provides global uniqueness
when processing a molecular collection. Call `clear_unique_molecule_cache()` to
forget previously generated structures. Clearing before each call gives
per-parent uniqueness while retaining the configured replacement collection.

The following controls are currently exposed:

- `set_ring_atom_smarts(smarts)` selects rings in the parent molecule.
- `set_unique_molecules_only(value)` controls duplicate suppression.
- `clear_unique_molecule_cache()` forgets products retained for duplicate suppression.
- `set_min_support_requirement(count)` sets the minimum precedent count.
- `set_max_formula_difference(value)` limits formula divergence.
- `set_remove_isotopes(value)` controls removal of attachment labels.
- `set_preserve_connectivity(value)` controls whether replacement records must
  preserve an observed substitution pattern; call it before reading rings.
- `read_replacement_rings(filename)` reads the replacement collection.
- `number_replacement_rings()` reports the retained collection size.
- `process(molecule)` generates and returns products.

For global uniqueness across a stream, configure and load the object once, then
reuse it without clearing:

```python
for molecule in molecules:
  for product in replacement.process(molecule):
    # Each structure is returned at most once across the stream.
    pass
```

To suppress duplicates independently for each parent:

```python
for molecule in molecules:
  replacement.clear_unique_molecule_cache()
  for product in replacement.process(molecule):
    # Duplicates are suppressed only among products of this parent.
    pass
```

The Python binding does not currently expose `ring_replacement_inexact`.

## Using `ring_replacement_inexact`

`ring_replacement_inexact` performs topology-changing scaffold hopping. It can,
for example, replace a fused `5a6a` aromatic system with a fused `6a6a` system,
or replace a fused bicyclic system with a single ring. The replacement need not
have the same number, sizes, or arrangement of rings as the parent scaffold.

The tool removes the selected ring system, retains its external substituents,
and enumerates ways of attaching those substituents to the replacement system.
This resolves the ambiguity introduced by changing topology, but can generate
many regioisomers and many more products than exact replacement.

### Quick start

The ring system to remove must be selected with `-s` or `-q`. For example:

```shell
ring_replacement_inexact \
  -R "${LILLYMOL_HOME}/data/ring_replacement/rings_6a6a.smi" \
  -s '[/IWfss2cr5r6]' -z i -e -n 5 -x 500 -I . -V -c \
  input.smi > variants.smi
```

This selects fused aromatic systems containing five- and six-membered rings and
tries fused `6a6a` replacements. `-e` retains same-ring attachment relationships
where possible, while `-x 500` limits combinatorial growth.

### Attachment placement and `-e`

By default, substituents are enumerated across allowed attachment positions on
the replacement scaffold. If several substituents were attached to the same
ring in the starting fused system, their relationship to one another may be lost
when the topology changes.

`-e` requires substituents that shared a ring in the starting scaffold to share
a ring in the product scaffold. It preserves this same-ring relationship, not an
exact ortho, meta, or para distance. Relationships involving different starting
rings are not constrained. Some transformations, such as reducing a fused
system to one ring, satisfy the same-ring requirement automatically.

Replacement records normally identify permitted attachment atoms. `-X any`
allows attachment to any available atom in the replacement ring system, further
increasing coverage and risk.

### Controlling product growth and risk

- `-n <count>` excludes replacement systems with fewer than `<count>` source
  examples.
- `-e` preserves same-ring attachment relationships from the parent.
- `-x <count>` approximately limits the number of products generated from each
  starting molecule. Enumeration checks the limit periodically, so it may be
  exceeded slightly.
- `-o` suppresses products with ortho substituent placement.
- `-V` discards products with invalid valence.
- `-Y <query>` requires complete products to match a query.
- `-N <query>` rejects complete products that match a query.

Because topology-changing replacements can produce a combinatorial number of
attachment arrangements, begin with a meaningful `-n` threshold, use `-e` when
same-ring relationships matter, and set `-x` during exploratory runs.

### Other options

- `-R <file>` reads a replacement collection and may be repeated. Use
  `-R F:<file>` when `<file>` contains a list of replacement-ring filenames.
- `-s <SMARTS>` and `-q <query>` identify the ring system to remove; at least one
  selection query is required.
- `-z i` ignores input molecules that do not match a selection query.
- `-p` writes the parent molecule before its products.
- `-I .` removes isotopic attachment labels from products. The argument is
  required; its value is currently not otherwise interpreted.
- `-c`, `-l`, and `-g ...` provide chirality removal, largest-fragment reduction,
  and chemical standardization.
- `-A ...`, `-E ...`, and `-i ...` provide standard aromaticity, element, and
  input-type controls.
- `-v` enables progress and summary information.

The `-F`, `-D`, `-d`, and `-f` options described for exact replacement are not
options of `ring_replacement_inexact`.


## Building replacement-ring collections

Most users can use the ChEMBL-derived collections distributed under
`${LILLYMOL_HOME}/data/ring_replacement`. Build a new collection when proprietary
chemistry, a different source corpus, alternative standardization, or atom-typed
attachment context is required.

Collection construction has two stages:

1. Run `ring_extraction` over each molecular corpus.
2. Optionally combine corresponding files from several corpora with
   `aggregate_replacement_rings`.

Extraction can be expensive for very large corpora but is performed only when a
collection is created or refreshed, not during normal ring replacement.

### Extracting rings with `ring_extraction`

A representative extraction is:

```shell
mkdir chembl_rings
cd chembl_rings

ring_extraction \
  -A 2 -S rings -R 7 -Z 3 -k -c -l -g all -r 100000 -v \
  /path/to/chembl.smi
```

This creates files named `rings_<class>.smi` in the current directory. Important
controls are:

- `-S <stem>` sets the required output stem; `-S rings` produces names beginning
  with `rings_`.
- `-R <size>` sets the largest individual ring size retained. The default is 7.
- `-Z <count>` sets the largest number of rings in a ring system. The default is
  3. Increase it when larger fused or spiro systems are important.
- `-k` also writes replacement records that do not constrain the observed
  substitution pattern. These records are used by `ring_replacement -d`.
- `-P <atype>` stores the atom type of each exocyclic atom attached to a ring.
  The same `-P` specification must be supplied to `ring_replacement`.
- `-a` changes aliphatic double bonds within rings to query bond type “any”.
- `-c`, `-l`, and `-g ...` remove chirality, retain the largest fragment, and
  standardize molecules before extraction.
- `-r <count>` reports progress periodically.
- `-X list=<file>` writes the names of all generated files to `<file>`.
- `-X noarom` generates aromaticity-independent queries, allowing aromatic and
  aliphatic ring forms to match one another.
- `-X sss` checks generated queries by searching the source molecule. This is
  primarily a validation option.

`-A 2` is a LillyMol aromaticity setting commonly used for these collections. It
allows rings with two pi electrons, notably some four-membered rings, to be
handled as aromatic.

Use consistent aromaticity, standardization, element, and atom-typing settings
when generating collections that will later be combined.

### Replacement-file naming

The class suffix describes the sizes and aromaticity of the component rings. On
case-sensitive filesystems, the default follows the SMARTS convention:

- lowercase `a` means aromatic;
- uppercase `A` means aliphatic.

For example:

| Filename | Contents |
| --- | --- |
| `rings_5a.smi` | five-membered aromatic rings |
| `rings_5A.smi` | five-membered aliphatic rings |
| `rings_5a6a.smi` | fused systems containing aromatic five- and six-membered rings |
| `rings_3A6a6A.smi` | systems containing aliphatic 3-, aromatic 6-, and aliphatic 6-membered rings |

Component ring descriptors are placed in a canonical order, primarily by ring
size. Spiro and fused systems use the same filename convention; the filename
describes the component rings, not their fusion type.

The default names are unsuitable for the commonly used case-insensitive macOS
filesystem because `rings_5a.smi` and `rings_5A.smi` cannot coexist. Generate
portable names by specifying:

```shell
ring_extraction -X Ar:Al -S rings ... input.smi
```

With this convention, `Ar` means aromatic and `Al` means aliphatic:

| Case-sensitive name | Case-insensitive name |
| --- | --- |
| `rings_5a.smi` | `rings_5Ar.smi` |
| `rings_5A.smi` | `rings_5Al.smi` |
| `rings_5a6a.smi` | `rings_5Ar6Ar.smi` |

The distributed data directory includes `to_mac.sh` and `to_linux.sh` for
converting an existing set of filenames. These scripts rename files; they do not
change their contents.

### Administrator wrapper

`contrib/bin/ring_replacement.py` is a small front-end intended for administrators
who maintain replacement-ring collections. It is not the primary user interface
for ring replacement. The wrapper lets an administrator choose a collection and
ring class by logical name and then constructs the corresponding `-R` arguments
for `ring_replacement`.

The `--ring` option uses canonical logical ring names regardless of the local
filesystem naming convention. Names are one to three component descriptors, each
consisting of a ring size followed by `a` for aromatic or `A` for aliphatic.
Components must be sorted by increasing ring size; when sizes are equal, `a`
precedes `A`. Examples include:

| Input | Meaning | Status |
| --- | --- | --- |
| `6a` | six-membered aromatic ring | canonical |
| `5a6a` | fused aromatic five- and six-membered rings | canonical |
| `6a6A` | aromatic and aliphatic six-membered rings | canonical |
| `6A6a` | same rings as `6a6A` | rejected; wrapper suggests `6a6A` |
| `6a5a` | same rings as `5a6a` | rejected; wrapper suggests `5a6a` |

On a case-sensitive filesystem, the wrapper looks for files such as
`<collection>_6a.smi` and `<collection>_6a6A.smi`. On a macOS-style collection
using `Ar`/`Al` filenames, the same user input is translated automatically; for
example `--ring 6a6A` resolves to `<collection>_6Ar6Al.smi`. The wrapper probes
for the common six-membered aromatic collection, `6a` or `6Ar`, to determine the
preferred filename convention, with a fallback based on the files present.

`--list-rings` reports canonical logical names, even if the underlying files use
`Ar`/`Al` names. This keeps administrator commands the same on Linux and macOS
collections.

### File format

Despite the historical `.smi` suffix, each file contains one
`RplRing::ReplacementRing` textproto per line, not an ordinary SMILES file. The
suffix remains for compatibility with existing workflows.

A record, formatted across several lines here for readability, resembles:

```textproto
smi: "[1CH]1=[1CH]C(=O)[1CH]=C[1NH]1"
smt: "[ax2r6D>2]1:[ax2r6D>2]:[ax2r6D3](:[A]):[ax2r6D>2]:[ax2r6D2]:[ax2r6D>2]:1"
id: "CHEMBL4552407.6a"
n: 186
conn: true
usmi: "O=c1[1cH]c[1nH][1cH][1cH]1"
```

The fields are:

- `smi`: the replacement ring structure.
- `smt`: the query used to locate a compatible parent ring. Its atom order is
  aligned with `smi`, enabling one-for-one replacement.
- `id`: the first source molecule retained as an exemplar.
- `label`: an optional descriptive class field; current extraction primarily records the class in the filename.
- `n`: the number of source molecules containing this ring.
- `conn`: whether the record constrains observed attachment connectivity.
- `usmi`: a canonical key used to combine equivalent records across corpora.
- `exo`: optional isotopically labelled atoms used to reconstruct exocyclic
  double bonds.

When `-k` is used, extraction also creates records without attachment-degree
constraints. Their `conn` field is false, and they are ignored by default exact
replacement but become available with `ring_replacement -d`.

### Spiro systems and exocyclic double bonds

`ring_extraction` processes spiro as well as conventionally fused ring systems.
The class name records the component ring sizes and aromaticity, while the record
itself contains the topology needed for replacement.

Exocyclic double bonds require additional information because including the
external atom directly in `smt` could prevent the query from matching a simpler
parent ring. Such atoms are stored separately in `exo`. Isotopes in `exo`
correspond to atom-map numbers in `smi`, allowing `ring_replacement` to restore
the double bonds when constructing a product.

### Extracting several source collections

Run unrelated corpora into separate directories so their counts and intermediate
files remain distinct:

```shell
mkdir corporate_rings chembl_rings

(cd corporate_rings && \
  ring_extraction -A 2 -S rings -k -c -l -g all -v \
  /path/to/corporate.smi)

(cd chembl_rings && \
  ring_extraction -A 2 -S rings -k -c -l -g all -v \
  /path/to/chembl.smi)
```

Use the same filename convention in each directory. If atom typing is required,
use the same `-P` specification for every corpus.

### Aggregating collections

`aggregate_replacement_rings` combines equivalent records by `usmi`, sums their
support counts, and retains the first exemplar encountered. Aggregate
corresponding ring classes, one output file at a time:

```shell
mkdir combined_rings

aggregate_replacement_rings \
  -S combined_rings/rings_6a.smi -v \
  corporate_rings/rings_6a.smi \
  chembl_rings/rings_6a.smi
```

Repeat this for each desired ring-system class. Do not combine different classes
such as `5a` and `6a` into the same output file. The resulting files can be passed
directly to `ring_replacement` with `-R`.
