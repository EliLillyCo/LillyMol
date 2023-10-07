# minor_changes

## TLDR

Want a potentially large number of new molecules that are derived
from a set of starting molecules:
```
minor_changes.sh -c file.smi > file.variants.smi
```

Probably generates more than you need. 500 random molecules from
ChEMBL took 30 minutes to generate 10.5M unique variants. There are
many ways in which these numbers can be reduced.

## Introduction

This tool is designed for making minor modifications to an existing molecule,
using simple transformations, as well as externally derived fragments. Many
attempts are made to have it generate plausible molecules, so the results
should mostly be reasonable, but this is definitely not guaranteed. Passing the
results through a synthetic feasibility assessment would generally be
desirable. Experience says that this may eliminate 80% of what is generated,
but this is obviously very dependent on the parameters selected.

## Complementary Tools
`ring_replacement` can be used for exploring different substitution patterns
in rings. By default, `minor_changes` does not swap atoms in an aromatic ring. Pipe
the output from `ring_replacement` to this tool in order to explore
further variants involving different ring atom arrangements.

Isostere replacement can be done with the reactions associated with the
`molecular_variants` tool - although some of those need to have the
reverse transformation also implemented. Again, a pipeline of generators
can be formed.

## Numbers
This tool can generate large numbers of molecules - largely driven
by the libraries of fragments specified, as noted previously. This
is despite the fact that the fragment libraries contained no more
than 4 atoms per molecule.

Again, it is recommended that the output be passed to a synthetic precedent
tool to eliminate unlikely atomic arrangements.

## Specifics.
The tool is complex, with a lot of choices made about how it behaves. For
that reason, the primary means of communicating with the tool is via a
`minor_changes_data::MinorChangesData` text proto, which is specified via the
`-C` option.

Here is one that I have found useful for testing.
```
# Allow removal of fragments containing as many a 3 atoms
remove_fragment: [1, 2, 3]
# Avoid excessive numbers of products for any starting molecule.
max_variants: 10000
# Running with atom typing is generally recommended
# This type contains the element and whether it is aromatic or not.
atype: "UST:AY"
# Do not use rare fragments. Helps cut numbers too.
# in this example, prefer setting library size
# fragment_support: 10

# Turn on all transformations

add_fragments: true
replace_terminal_fragments: true
single_to_double_bond: true
double_to_single_bond: true
unspiro: true
make_three_membered_rings: true
change_carbon_to_nitrogen: true
change_carbon_to_oxygen: true
change_nitrogen_to_carbon: true
insert_ch2: true
remove_ch2: true
destroy_aromatic_rings: true
destroy_aromatic_ring_systems: true
swap_adjacent_atoms: true
swap_adjacent_aromatic_atoms: true
insert_fragments: true
replace_inner_fragments: true

# These substantially cut the run time.
max_fragment_lib_size: 100
max_bivalent_fragment_lib_size: 200
```
The proto definition is [GitHub](https://github.com/EliLillyCo/LillyMolPrivate/blob/a6b84fa94d451438a6a16166c9eaf4b5f4c76d53/src/Molecule_Tools/minor_changes.proto#L1)

The instance above takes 52 seconds to generated 680k variants from 500 input
molecules.

## Atom Typing.
It is recommended that the tool be used with atom typing enabled. As
seen in the proto above, I have been using `UST:AY` which classifies
atoms by their atomic number and aromaticity only. Other atom types
would clearly be possible.

The reason to use atom types in the library is that ensures that the
fragment is joined to an atom similar to the atomic context from
which it was extracted. For example a phenolic Oxygen atom would
record an atom type of `aromatic carbon`, and so that particular
oxygen would only be attached to aromatic carbons. There will be
other single atom oxygen fragments that record that they were 
attached to an aliphatic carbon... This should contribute to better
fragment generation.

### Inputs
There are two kinds of external structure files that can be provided, which
provide fragments of various kinds.

* Monovalent Substituents/Replacements
* Bivalent  Substituents/Replacements

Monovalent fragments are joined to atoms with complementary atom types
that have an available Hydrogen atom; `add_fragments` directive.

Monovalent fragments can also replace an existing substituent of size
specified with the `remove_fragment` directive. The
existing terminal group is removed, and the fragment joined, if the
atom types match, to the join point.  `replace_terminal_fragments` directive.

Bivalent fragments can have either 1 or 2 attachment points. If there is
a single attachment point, it is inserted between two atoms, and each of
those two atoms bond to the same atom in the fragment. If there are two
attachment points, a bond in the parent molecle is selected and removed.
The bivalent fragment is inserted by bonding (in both directions) to the
remaining atoms in the parent.

You can limit the number of library items used via the `max_fragment_lib_size`
and `max_bivalent_fragment_lib_size` values in the proto. This can
substantially cut run times, but at the expense of not producing what might
be quite plausible variants.

## Generating Library files.

Skip if you just want to use the tool. But if you find you are getting
too many molecules generated, this is likely the cause.

The tool reads modified `dicer_data::DicerFragment` protos. The modification
is that the smiles is written as the first token on the line, folowed by the
proto. This has the desirable effect of making the file readily processed by
all standard tools, at the cost of requiring special processing in order
to read the proto.

It is possible that `dicer` could generate both the monovalent and the
bivalent substitunts, here I show using `get_substituents` and `dicer`.

### Monovalent Fragments
Use `get_substituents` to extract fragments from known collections. In this
invocation, we are looking for fragments with a max of 3 atoms.
```
get_substituents.sh -s '*' -M 3 -n -P UST:AY -S S3.txt -v chembl.smi proprietary.smi ...
```
The smarts `*` means any atom, so we look for all contexts for substituents.

The resulting file S3.txt might look like
```
O=[3038NH] iso: ATYPE smi: "O=[3038NH]" par: "CHEMBL9804" nat: 2 n: 567 
ON=[6007CH2] iso: ATYPE smi: "ON=[6007CH2]" par: "CHEMBL101180" nat: 3 n: 105 
Br[21001CH]=C iso: ATYPE smi: "Br[21001CH]=C" par: "CHEMBL17255" nat: 3 n: 2 
[3001NH3] iso: ATYPE smi: "[3001NH3]" par: "CHEMBL153534" nat: 1 n: 316390 
[24001CH]#CC iso: ATYPE smi: "[24001CH]#CC" par: "CHEMBL3350642" nat: 3 n: 2 
F[9001CH2]N iso: ATYPE smi: "F[9001CH2]N" par: "CHEMBL3897128" nat: 3 n: 1 
```
Again, note the extra token in column 1. The isotopes are the atom types
of the atoms to which these fragments used to be attached. These numbers of
course make no sense on their own. If curios, use `fileconv` to get molecules
labelled by atom type
```
fileconv.sh -I atype -Y atype=UST:AY -S - file.smi
```
which specifies an atom typing of `UST:AY` and that the isotopic label is
the atom type. This will show that isotope 3001 corresponds to an
aliphatic Carbon atom. Normally you will not need to pay any attention
to these numbers, just use the atom type specification.

If I were to make the atom typing more precise, the next step would be
to differentiate saturated from unsaturated atoms, `UST:ABY`.

`get_substituents` is way faster than `dicer`, processing both ChEMBL and
the corporate collection in about 4 minutes, whereas `dicer` takes close to
an hour to do the same molecules.

Do any further filtering on this file to eliminate fragments you do not want.

See the script `get_substituents_to_fragments.sh` in the data directory for an example
of such filtering.

### Bivalent Fragments
`dicer` can be used to generate bivalent (and monovalent) fragments. A typical
usage might be
```
dicer.sh -c -P UST:AY -B nooutput -B nbamide -B brcb -X 500 -M 4 -I atype -z i -B smiles_proto -B fragstat=fragstat_atype.txt chembl.smi ...
```
The output from this might look like
```
[9001NH2]CC[3038NH2] smi: "[9001NH2]CC[3038NH2]" par: "CHEMBL22077" nat: 4 n: 222 
[3038SH][9001NH][3038CH3] smi: "[3038SH][9001NH][3038CH3]" par: "CHEMBL53051" nat: 3 n: 10 
[3038CH3][6007CH]1[9001CH2]C1 smi: "[3038CH3][6007CH]1[9001CH2]C1" par: "CHEMBL1370756" nat: 4 n: 5 
F[3001CH2][9001CH3] smi: "F[3001CH2][9001CH3]" par: "CHEMBL1822932" nat: 3 n: 6 
[6007CH2]=[6007CH]C smi: "[6007CH2]=[6007CH]C" par: "CHEMBL75824" nat: 3 n: 120 
[3038CH3][6007CH2][3038CH3] smi: "[3038CH3][6007CH2][3038CH3]" par: "CHEMBL11390" nat: 3 n: 92 

```
Again, this almost certainly needs to be filtered to remove fragments
not wanted in enumerated variants.

See the script `dicer_frags_to_bivalent.sh` in the data directory for an example
of such filtering.

There will also be fragments that do not have a free hydrogen atom on the
labelled atom. This is usually because of Sulphur atoms, which have a variety
of valence states, for example `O=[9001S]=O`. The filtering scripts above
will get rid of these.

Note too that the protos also include the number of exemplar molecules. At run
time the `-u` option can be used to impose a support requirement on the
libraries being read, which can significantly reduce the numbers. You will
see messages about failure to build a fragment for those fragments falling
below the support requirement.

## Details
Internally the tool defines several kinds of molecular transformations
which are applied exhaustively to all plausible atoms. For example the
transform Carbon to Nitrogen transformation can only be applied
to Carbon atoms. But in order to avoid creating implausible molecules, it
will not change a Carbon to
a Nitrogen if the Carbon atom already has a heteroatom attached.

The transformations are ordered so that those which can generate the
largest numbers of variants are done last, so that if you impose
a limit on the number of variants per starting molecule, the most
prolific variants are done last.

Note that even if a max number of variants is specified, you may get more
than that number produced, because the tool only checks periodically on
the number if items generated.
