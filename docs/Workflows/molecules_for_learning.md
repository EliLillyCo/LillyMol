# Generating Molecules for Pretraining

This HOWTO describes LillyMol workflows for generating virtual molecules from a
starting set of known molecules. The examples are aimed at data augmentation and
pretraining workflows, where the goal is not to enumerate every possible analogue
but to generate a large, chemically plausible population around molecules of
interest.

Chemprop's Chemeleon pretraining work has shown that pretraining can be useful.
The question addressed here is how LillyMol tools can generate additional
molecules that might be useful for similar learning tasks.

The workflows below assume that you already have a starting SMILES file,
`start.smi`, containing molecules worth augmenting.

## Choosing a Generator

Several LillyMol tools can generate new molecules from a starting set. They use
very different strategies, and it is usually worth trying more than one method.

| Tool | Strategy | Typical scale |
| ---- | -------- | ------------- |
| [medchem_wizard](/docs/Molecule_Tools/medchemwizard.md) | reaction-driven transformations | moderate |
| [minor_changes](/docs/Molecule_Tools/minor_changes.md) | small local edits to one molecule | moderate |
| [sidechain_switcheroo](/docs/Molecule_Tools/sidechain_switcheroo.md) | swaps sidechains found within the input set | high, input-size dependent |
| [random_molecular_permutations](/docs/Molecule_Tools/random_molecular_permutations.md) | larger random edits from fragment rules | tunable |
| [safe_generate](/docs/Molecule_Tools/SAFE.md) | text-based SAFE fragment replacement and breeding | high, library dependent |
| [ring_replacement](/docs/Molecule_Tools/ring_replacement.md) | exact ring replacement from a ring library | moderate |
| [ring_replacement_inexact](/docs/Molecule_Tools/ring_replacement.md) | broader ring replacement from a ring library | moderate to high |
| [substituent_identification](/docs/Molecule_Tools/substituent_identification.md) | local matched-pair-like substituent replacement | library dependent |

Tools that make intermolecular combinations, such as `sidechain_switcheroo` and
`safe_generate`, can grow rapidly with the diversity and size of the input set.
For those tools, splitting the input can change the products generated. For tools
that operate primarily on each molecule independently, splitting large inputs is
usually a good way to manage long runs.

The tools do not all produce identical output conventions. Product names may
contain different annotations, duplicate suppression varies by tool, and some
outputs may contain residual isotopic labels. Treat the raw output as an
intermediate file and use a common cleanup step before combining products.

Many generators will briefly form products with impossible valence states while
searching. Most such products are discarded, but warnings can still be noisy. For
large runs, redirect standard error to a log file.

```bash
command ... input.smi > products.smi 2> products.log
```

All tools have many additional options. The commands below are starting points;
consult the linked tool documentation before large production runs.

## Immediate Generators

These tools are the easiest to run because they do not require a precomputed
fragment or ring database.

### Medchem Wizard

`medchem_wizard` applies reaction-like medicinal chemistry transformations.

```bash
medchem_wizard.sh -U all -m 20 -M 50 -c -V . -y -l -v start.smi > products.smi
```

Adjust the lower and upper atom-count thresholds with `-m` and `-M`. The
`-U all` option suppresses duplicate products across all starting molecules;
this is useful, but memory use increases with the number of products retained.

### Minor Changes

`minor_changes` makes small, chemically conservative changes to each input
molecule.

```bash
minor_changes.sh -x start.smi > products.smi
```

This tool currently does not discard duplicate products, so post-process with
[unique_molecules](/docs/Molecule_Tools/unique_molecules.md) or the common
cleanup workflow below.

```bash
minor_changes.sh -c -l -x start.smi | unique_molecules -g all -S - - > products.smi
```

### Random Molecular Permutations

`random_molecular_permutations` makes more substantial changes than
`minor_changes`. It uses installed fragment libraries, and additional libraries
can be supplied when needed.

```bash
random_molecular_permutations.sh -r 5 -R 7 -c 15 -C 50 -y 1 -Y 6 -G +1 -G -1 start.smi > products.smi
```

This example constrains the generated molecules by ring count, atom count, and
other molecular properties. Tighten these limits if the output becomes too broad
for the intended learning task.

## Collection-Dependent Generators

These tools use relationships among molecules or fragments. They can generate
many more products, but the product set depends strongly on the input collection
or on a prebuilt library.

### Sidechain Switcheroo

`sidechain_switcheroo` identifies sidechains and swaps them among molecules in
the input set.

```bash
sidechain_switcheroo -s '[R]!@-*' -z i -x 2 -I -V -c -v start.smi > products.smi
```

This can generate very large numbers of molecules from large, diverse input
sets. Because the fundamental operation is exchanging sidechains within the
input, splitting the input changes the accessible sidechain combinations and is
not advised unless that is intentional.

### SAFE Generate

`safe_generate` works on SAFE SMILES, where isotopically marked connection
points make fragment replacement a text operation. First convert the starting
molecules to SAFE form and write a fragment-library summary.

```bash
mol2SAFE -I 1 -S start.lib.textproto start.smi > start.safe.smi
```

It is possible to run only on the starting set, but output is usually more
interesting when external fragments and linkers are also supplied. LillyMol ships
with dicer fragment summaries that can be converted into a SAFE fragment library.

```bash
dicer_fragments_collate -p 10 -nosmi -v \
        ${LILLYMOL_HOME}/data/dicer_fragments/FRAG_1_[1-5].textproto > frag.textproto

fileconv -S - frag.textproto | mol2SAFE -I 1 -S frag.lib.textproto - > frag.safe.smi
```

This gathers single-attachment-point ChEMBL fragments up to five heavy atoms and
writes `frag.lib.textproto`, a library that `safe_generate` can use.

```bash
safe_generate -n 1000 -b 1000 -e 0 -C safe_generate.textproto \
  -L frag.lib.textproto -L start.lib.textproto start.safe.smi > products.smi
```

In this example, `-n 1000` requests up to 1000 products from each starting
molecule by replacing SAFE fragments with library fragments. `-b 1000` requests
up to 1000 products by breeding fragments between input molecules. The `-e 0`
setting disables exhaustive library replacement; increase it if you want that
mode.

### Substituent Identification

`substituent_identification` is a local matched-pair-like generator. During a
database build, it fragments known molecules and stores substituents keyed by
the circular-fingerprint context in which each substituent occurred. During
lookup, a starting molecule is fragmented in the same way and compatible
substituents from the database are proposed as replacements.

Build the database from a large collection of real molecules, such as ChEMBL or
an internal collection.

```bash
substituent_identification -B -d chembl.frags.bdb -R 5 -w 10 -M 12 -v \
  -Y dbproto -Y rpt=10000 chembl.smi
```

This may take a while. The `-Y rpt=10000` directive reports progress every
10,000 molecules. The database size depends strongly on the maximum radius
specified with `-R` and the maximum substituent size specified with `-M`.

Once the fragment database has been built, use it to replace substituents in the
starting set. This example breaks a bond between a ring atom and an exocyclic
atom and replaces the removed fragment with database fragments seen in similar
contexts.

```bash
substituent_identification -d chembl.frags.bdb -s '[r]-!@*' -k -I \
  -Y maxgen=100 start.smi > products.smi
```

Use more specific SMARTS when you want to replace only particular substituents.
The `-I` option removes isotopic labels from products before writing. The
`-Y maxgen=100` limit is a useful guard for broad queries and large databases.

## Ring Generators

These tools replace ring systems using curated ring libraries. They are often a
good way to generate realistic scaffold variants, but their products may look
farther away by fingerprint similarity than the chemical change suggests.

### Ring Replacement

`ring_replacement` performs precise replacement of rings matched by a query. A
very generic ring query plus a broad replacement-ring library can generate many
plausible variants.

The file [common_rings](/data/ring_replacement/common_rings.txt) contains
replacement rings that can be used with both `ring_replacement` and
`ring_replacement_inexact`.

```bash
ring_replacement -R F:${LILLYMOL_HOME}/data/ring_replacement/common_rings.txt -s '[R]' \
        -z i -u -d -I . start.smi > products.smi
```

### Ring Replacement Inexact

`ring_replacement_inexact` allows broader replacement than exact ring
replacement.

```bash
ring_replacement_inexact -R F:${LILLYMOL_HOME}/data/ring_replacement/common_rings.txt -s '[R]' \
        -z i -I . start.smi > products.smi
```

If you plan to run `ring_replacement_inexact`, running exact `ring_replacement`
may be unnecessary because exact replacements should generally be a subset of
the inexact output.

### Others
[trxn](/docs/Molecule_Tools/trxn.md) performs reactions. For a given set of
starting molecules there may be project specific transformations that could
be applied, generating variants of the starting set.

## Common Cleanup

Generated products from different tools can be identical. They can also contain
residual isotopic labels, multiple fragments, or structures that should be
rejected by normal LillyMol standardisation and valence checks.

Run each raw product file through `fileconv` before combining product streams.

```bash
fileconv -V -g all -I change -v -f lod -S products_ok products.smi
```

This removes structures with valence errors, applies chemical standardisation,
keeps the largest fragment, and changes isotopic atoms to non-isotopic form. You
may also want atom-count limits with `-c` and `-C`, or ring limits with `-r` and
`-R`.

If you want to retain the generator identity, add that label before combining
all product files.

These virtual molecules were not generated with a general reactivity or
undesirable-functionality filter, so this is also a good stage to apply Medchem
Rules or project-specific filters.

After cleanup, combine all product files and deduplicate them. For large files,
use `unique_molecules_parallel.sh`.

```bash
unique_molecules_parallel.sh -thr 16 -S unique generated.smi
```

Use a thread count appropriate for the machine. Depending on how many molecules
were generated, this may be the longest step in the workflow.

## Proximity to Starting Molecules

Some workflows require generated molecules to remain close to the starting set.
This is project-dependent. For example, exact ring replacements may be chemically
conservative even when a fingerprint Tanimoto distance changes substantially.
`ring_replacement` also has options that limit molecular formula changes during
replacement.

Use `gfp_distance_filter` to remove products that are too far from every
starting molecule.

```bash
gfp_make.sh starting_molecules.smi > starting_molecules.gfp
gfp_make.sh products.smi | gfp_distance_filter -T 0.25 -N 1 -f - > okdistance.smi
```

This keeps products within 0.25 distance of at least one starting molecule. The
right threshold depends on the objective of the enumeration and should be tuned
empirically.

## Summary

The generators discussed here take different approaches to molecule generation.
Some operate on each molecule independently, some depend on the whole input
collection, and some require prebuilt fragment or ring libraries. Expect overlap
between methods, and plan on a common cleanup and deduplication pass before using
the generated molecules for training.
