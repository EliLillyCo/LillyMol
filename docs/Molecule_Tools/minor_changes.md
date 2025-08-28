# minor_changes

## Objective

This tool uses rule based transformations to make small changes to starting
molecules. The most common usage is to help fill out an SAR around
a set of active molecules. The transformations applied might
include adding small fragments, removing CH2 groups, changing
a substitution pattern, etc... Since it applies all available
transformations exhaustively, potentially large numbers of molecules
can be formed.

Taking the default, which turns on all known transformations that
do not depend on external fragments, and which operates on all
possible sites in the starting molecule
```
minor_changes file.smi > file.variants.smi
```
generates 31k molecules from a starting set of 1000 in 1.2 seconds.

The rules applied can be controlled via a configuration file, which will
also enable specification of fragments to be added, and reactions to be
performed. A shell script [minor_changes.sh](../../contrib/bin/minor_changes.sh)
invokes the executable with all options turned on, as well as making
available some reaction based transformations, as well as a small
fragment library (60 functional groups, max 3 atoms). Given
1000 random starting molecules the tool generates 337k new molecules in about
15 seconds using that configuration file.

Obviously if a larger fragment library is used, those number can increase
dramatically - see below for how to generate fragment libraries.

The reactions specified with the script include a number of common
isostere type transformations. More could be added - definitely open
to ideas... Reactions in the [medchem wizard](/data/MedchemWizard) tool
could also be used here.

The tool struggles with chirality and needs fixes to make that robust.
For now, the safest thing to do is remove chirality from your input
molecules via the -c option.

## Details

The tool is designed for making minor modifications to an existing molecule,
using simple transformations, as well as externally derived fragments.  Many
attempts are made to have it generate plausible molecules, so the results should
mostly be reasonable, but this is definitely not guaranteed.  Passing the
results through a synthetic feasibility assessment would generally be desirable. 
Experience says that this may eliminate a significant fraction of what is
generated, but this is obviously very dependent on the parameters selected and
the starting molecules.

## Complementary Tools
`ring_replacement` can be used for exploring different substitution patterns
in rings. By default, `minor_changes` does not swap atoms in an aromatic ring. Pipe
the output from `ring_replacement` to this tool in order to explore
further variants involving different ring atom arrangements.

[substituent_identification](substituent_identification.md) can perform 
sidechain replacements of varying precision.

[random_molecular_permuations](random_molecular_permuations.md) also contains
various transformations, possibly not as well developed as in minor_changes,
but which can often generate interesting modifications of a starting structure.

[smiles_mutation](smiles_mutation.md) performs random string based transformations
on smiles strings, attempting molecule interpretation on the resulting strings.
Most strings generated are invalid.

[SAFE](SAFE.md) provides a very flexible and powerful means of generating
molecular variants. Although a newer tool, developed after the publication
of SAFE, it is probably the most 'interesting' de-novo generator in LillyMol.


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
# This must be the same as what was used to build the library fragment files.
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
remove_fused_aromatics: true

# These substantially cut the run time.
max_fragment_lib_size: 100
max_bivalent_fragment_lib_size: 200
```
The proto definition is [GitHub](../../src/Molecule_Tools/minor_changes.proto)

## Atom Typing.
If library fragments are being added, it is recommended that the tool be used
with atom typing enabled.  As seen in the proto above, I have been using
`UST:AY` which classifies atoms by their atomic number and aromaticity only. 
Other atom types would clearly be possible.

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
attachment points, a bond in the parent molecule is selected and removed.
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
bivalent substitunts, here I show using `get_substituents` and `dicer`,
although it should be noted that if you need to run `dicer` in order to 
get bivalent fragments, you get the substituents essentially for free
during that run.

If you just care about between ring linkers, you might be able to use 
get_linkers in order to generate bivalent linkers.

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
course make no sense on their own. If curious, use `fileconv` to get molecules
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

More information about the transformations are in the [proto](../../src/Molecule_Tools/minor_changes.proto)
