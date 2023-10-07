# Ring Replacement

## Motivation
Given a molecule of interest, frequently there is interest in variants of
that molecule. While similarity searching can identify existing molecules
that might be similar to a molecule, that does not answer the question
of what plausible variants might exist?

`ring_replacement` takes a starting molecule, a database of known rings,
and generates molecules that have the existing rings replaced by rings
from the database. The way the database is constructed, ring sizes and
substitution patterns are preserved. So if the starting molecule
is a para substituted benzene ring, the results will also be para
substituted 6 membered aromatic rings.

For example of the starting molecule is
```
Clc1ccc(F)cc1 p-benzene
```
requesting replacement of the 6 membered aromatic
```
ring_replacement -R rings_6a.smi p-benzene.smi
```
generates 185 new molecules. The first 10 (from Chembl) of those are
```
Cl[1c]1cc[1c](F)cc1 p-benzene %% CHEMBL503634.6a 1496512
Cl[1c]1c[n][1c](F)cc1 p-benzene %% CHEMBL156037.6a 83606
Cl[1c]1[n]c[1c](F)cc1 p-benzene %% CHEMBL156037.6a 83606
Cl[1c]1[n]c[1c](F)c[n]1 p-benzene %% CHEMBL1171471.6a 19926
Cl[1c]1c[n][1c](F)[n]c1 p-benzene %% CHEMBL1171471.6a 19926
Cl[1c]1[n][n][1c](F)cc1 p-benzene %% CHEMBL600052.6a 10360
Cl[1c]1[n]c[1c](F)[n]c1 p-benzene %% CHEMBL268339.6a 7938
Cl[1c]1cc[1n+](F)cc1 p-benzene %% CHEMBL505408.6a 6908
Cl[1n+]1cc[1c](F)cc1 p-benzene %% CHEMBL505408.6a 6908
```
The output consists of the smiles of the new molecule.  The name of
the starting molecule.  Then follows the name of an exemplar molecule,
that contains an example of the ring that has been
inserted.  The `6a` suffix indicates that the replacement is a six
membered aromatic.  The last token is the number of these rings found
in the knowledge base.  The output is sorted by occurrence.

Note that due to symmetry, many of the replacement rings are used twice.

Clearly the larger the number of examples of a ring in the
knowledge base, the higher the probability that this ring
might be synthietically feasible. Experience tells us that rings 
with low numbers of exemplars should be treated with caution.
While some may represent real molecules that have been made,
many seem more likely to have been drawing errors, or aspirations.
For example:
```
Fc1[c+]c[n+](Cl)cc1 pbenzene %% PBCHM59575064.6a 2
Fc1ccc(Cl)[p][p]1 pbenzene %% PBCHM159136080.6a 2
Fc1c[nH][n](Cl)so1 pbenzene %% PBCHM144506470.6a 2
Fc1[c-][o+]c(Cl)cc1 pbenzene %% PBCHM139100082.6a 2
```
are rings that have only two examples in the knowledge base.

## Usage
In order to make ring replacement work, two procedures are required.
1. Scan one or more collections of known molecules and extract the rings.
2. Use those knowledge bases to perform replacements.

# ring_extraction
ring_extraction is the tool that gathers data on the existing collections.
The following options are recognised.

```
Extracts rings and ring systems creating ReplacementRing protos that can be used by ring_replacement
 -S <stem>         create ring data files with <stem>
 -R <rsize>        max ring size to process (def 7)
 -Z <size>         max ring system size to process (def 3)
 -k                also generate smarts with connectivity not specified
 -a                transform within ring aliphatic double bonds to type any\n";
 -P <atype>        label atoms by atom type of exocyclic attached atom
 -X ...            more options
 -c                remove chirality
 -g ...            chemical standardisation
 -l                strip to largest fragment
 -v                verbose output
```
A common usage might be
```
ring_extraction -A 2 -S collection -k -c -l -g all -v collection.smi
```
for several different collections (corporate, chembl, pubchem, ...). Might be a good
idea to do each one in a separate directory.

The `-A 2` option combination tells the tool to consider rings that have
only 2 pi electrons to be aromatic. This is usually 4 membered rings, that
can usefully be considered aromatic.

```
mkdir corporate
cd corporate
ring_extraction -A 2 -S rings -k -c -l -g all -v /path/to/corporate.smi
cd ..
mkdir chembl
cd chembl
ring_extraction -A 2 -S rings -k -c -l -g all -v /path/to/chembl.smi
cd ..
mkdir ...
```


## Naming convention.
We adopt an idea from smarts in that aromatic atoms are represented
as lowercase and aliphatic as uppercase. So the file containing data
on 6 membered aromatic rings will be `rings_6a.smi` and the file containing
data on 6 membered aliphatic rings will be `rings_6A.smi`.

In the case of fused rings, there will be multiple ring sizes and aromaticity,
so a benzimidazole ring would be found in `rings_5a6a.smi`. Names are always
sorted from smallest to largest ring. The file `3A6a6A.smi` contains
fused (or possibly spiro fused) systems consisting of a
. 3 membered aliphatic
. 6 membered aromatic
. 6 membered aliphatic

Such a dataitem might be
```
smi: "[1CH]1=C2CCC3(CC3)OC2=N[1CH]=N1"
     smt: "[ax2r6D3]1:[ax3r6r6D3]2:[ax2r6D2]:[ax2r6D2]:[ax4r3r6D4]3(:[ax2r3D2]:[ax2r3D2]:3):[ax2r6D2]:[ax3r6r6D3]:2:[ax2r6D2]:[ax2r6D3]:[ax2r6D2]:1"
     id: "SCHEMBL2706866" n: 126 conn: true usmi: "O1C2(CCc3c1[n][1cH][n][1cH]3)CC2"
```
which includes a spiro fused 3 membered ring.

Note that the default setting of imposing a maximum ring system size of 3
will exclude a common motif of four rings, with a spiro fusion in the
middle. Adjust the -Z option as needed.

The letters 'a' and 'A' can be changed via the `-X` option
to `ring_extraction`. This may be important on file systems where file names
are case insensitive, where the files 'rings_5a.smi' and 'rings_5A.smi' cannot
coexist. See the '-X' option.

On a fairly old computer (2014), running ring extraction on 2.2M Chembl molecules
takes 5 minutes. Newer hardware will see that done in under 2 minutes.

## Details - TLDR, for geeks only.

In each directory there will be a number of files of the form `rings_*.smi`.
The naming describes what kind of rings are in that file. For
example `rings_6a.smi` contains the 6 membered rings from the collection. That
file might contain entries that look like
```
smi: "[1CH]1=[1CH]C(=O)[1CH]=C[1NH]1" smt: "[ax2r6D>2]1:[ax2r6D>2]:[ax2r6D3](:[A]):[ax2r6D>2]:[ax2r6D2]:[ax2r6D>2]:1" 
     id: "CHEMBL4552407.6a" n: 186 conn: true usmi: "O=c1[1cH]c[1nH][1cH][1cH]1" 
```
where we have split what will be a single line in that file.

The magic of ring replacement is the first two tokens, 'smi' and 'smt'. These
are aligned with each other. So if 'smt' is used to identify the atoms in
a ring system, the atoms from 'smi' can replace those atoms one-per-one and
get a replacement ring that preserves the connectivity pattern.

The 'id' field will be the first molecule in the input that exemplified this
particular ring. The 'n' field contains the number of molecules in the
collection that exemplify this ring. The final 'usmi' field is a unique
smiles, whose only purpose is to enable merging data generated across multiple
collections.

Note that if the `-k` option was used you will also have entries like
```
smi: "O=N1=CC=NN=C1" smt: "[A]=[ax2r6]1:[ax2r6]:[ax2r6]:[ax2r6]:[ax2r6]:[ax2r6]:1"
     id: "CHEMBL1457335.6a" n: 12 usmi: "O=[n]1cc[n][n]c1" 
```
where the unique smiles does not contain isotopic labels. Notice too that
the query atoms of the second query do not contain 'D' directives and
so can match any substitution pattern. Note too that the 'conn' (connections) 
attribute is not set.

### Spiro Fused Rungs
Unlike other LillyMol tools, this one generates rings that span sipro
fusions. Normally this will be seamless, and if you specify a sprio fused
system to be replaced, it will be replaced. But spiro rings are not
called out specifically in how the rings are stored.

### Exocyclic Double Bonds
Earlier versions of this tool did not handle exocyclic double bonds. With this
new version, for example
'O=C1NC=CC=C1F CHEMBL4558322' can replace the ring in 'N(C)(C)C1=CC=NC=C1 CHEMBL3561645'
to generate 'N(C)(C)[1C]1=CC=CNC1=O CHEMBL3561645 %% CHEMBL4558322.6a'. This new
functionality may introduce some undesirable changes as well, but is generally
desirable.

The proto data for these rings with exocyclic double bonds will be somewhat
different, possibly looking like
```
smi: "[CH2:70]1NC=CC=[1CH]1" smt: "[ax2r6D2]1:[ax2r6D2]:[ax2r6D2]:[ax2r6D2]:[ax2r6D2]:[ax2r6D>2]:1"
        id: "CHEMBL4558322.6a" n: 1 conn: true exo: "[70O]" usmi: "O=c1[nH]ccc[1cH]1"
```
The exocyclic atoms are specified as a separate molecule with an isotopic label near
70. Somewhere within the replacement ring 'smi:', there will be an atom map number
corresponding to this isotopic label. That is where the doubly bonded exocyclic
atom will be joined. The corresponding atom within the smarts does still have 'D2'
so we can not yet do fully flexibile ring replacement. Getting there...

This distribution does include rings extracted from a recent version of Chembl. See
[data](/contrib/data/ring_replacement) where you will find some 'hidden' files.

This version also expands coverage to ring systems of size 3 by default, but
currently there is really no limit to what can now be processed - see the '-Z'
option.

## Aggregating Across Collections
There is a script that can be used to aggregate data from multiple collections
into a single collection.

```
ruby aggregate_rings.rb corporate/rings chembl/rings ...
```
will aggregate all the data in the individual collections and create a unified set
of protos in the current directory. The example used will be the first example
encountered during the collation process. The script is quite dumb
and mostly treats the protos as text, but it works well enough and
is fast enough.

## Ring Replacement
Once a set of replacement rings has been assembled, those can be used to
perform ring replacement on molecules with existing ring/ring systems.
 -R <fname>    file of labelled rings created by ring_extraction
 -s <smarts>   only replace rings matched by <smarts>
 -q <query>    only replace rings matched by <query>
 -u            unique molecules only
 -a            allow change of aromaticity on ring replacement
 -p            write parent molecule
 -n <ex>       only process replacement rings with <ex> or more examples
 -w            sort output by precedent count
 -d            do NOT preserve substitution patterns. Replacement rings may
               not have had the same substition pattern as they do here
 -Y <query>    product molecules MUST     match a query in <query>
 -N <query>    product molecules must NOT match any queries in <query>
 -D <query>    discard any replacement ring that matches <query>
 -I .          remove isotopes from product molecules
 -I <n>        change all existing isotopes to <n> (useful if atom types used)
 -B <fname>    write molecules not transformed to <fname>
 -X ...        miscellaneous options, enter '-X help' for info\n";
 -c            remove chirality
 -l            strip to largest fragment
 -v            verbose output
```

So if you wanted to replace a 6 membered aromatic ring in some starting molecules
that might be done via
```
ring_replacement -u -v -R RINGS_6a.smi -s '[/IWfss1cr6]' ... input.smi
```
which for 2k input molecules, generates 22k new structures in 18 seconds. Only
1450 molecules had an isolated six membered aromatic, so for those molecules
that could be changed, we generated an average of 15 new structures. Note that
the `-u` option significantly reduced the number of molecules generated by
removing duplicates, so the actual number of new variants per starting ring
is larger.

### Which ring
If specified, only rings matched by the -s or -q options are considered for
replacement. Once the query match is made, any ring that contains a matched
atom is considered for replacement.

### Filtering
The -Y and -N filters apply filters to the generated molecules. The same
effect could be achieved by passing the results to `tsubstructure`, but
this is more efficient.

### Precedent
There is no 'right' answer to what value to use. If the -n option is omitted
all replacement rings are considered, regardless of precedent. This might be
what you want. Again, rings with very few examples might be risky in terms
of whether the proposed new molecule can actually be made.

## Atom types.
You may wish to preserve the context of the ring. For example if the starting
molecule has an OH group, you might only want to consider replacements where
the donor/replacement ring also had an OH group at that position too. In that
case we need to build the exemplar collection with atom typing and then use that same
atom typing when doing replacements.

In this case, we use an atom typing that consists of just number of connections and
atomic number, although generally I think there might be better choices.

```
ring_extraction -k -c -A 2 -S RINGS -P UST:CY -l -v chembl.smi
```
```
ring_replacement -A 2 -R RINGS_5a6a.smi -P UST:CY flubendazole.smi
```
Things will silently fail, no matches, unless the atom types used during
database building and querying are the same. Clearly it would be desirable
to make this more robust.

A useful atom typing might be `UST:ABCHY`. This consists of five atomic properties,
and as such, is likely over-specified. The comonents are

* A aromatic
* B unsaturation (excluding aromatic)
* C number of connections
* H number of hydrogens
* Y compressed atomic number - heavy halogens equivalent.

Using a very precise atom typing like this should raise the probability of your
newly generated molecules being synthetically feasible, but at the cost of generating
potentially many fewer new molecules.

Atom typing is described in [atom typing](/docs/Molecule_Lib/atom_typing.md).
