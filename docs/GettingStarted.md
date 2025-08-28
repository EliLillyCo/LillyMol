# LillyMol
LillyMol command line tools support a wide variety of commonly performed
cheminformatics related tasks.

The fundamental unit of work will be one or more files containing chemical
structures.

If you start with structures in a .sdf file, you probably need to start with
fileconv to convert that to smiles, using options to transfer information from
the SDF tags to the smiles. For example a download from Enamine stores the
ID of the molecule in the 'ID' tag.
```
>  <ID>
EN300-07843

```
We can transfer this to the molecule name with
```
fileconv -i sdf -i SDFID:ID -v .... file.sdf
```
which generates output like
```
CSC(C)CC=O ID:EN300-7099459
```
If we don't want the 'ID' token
```
fileconv -i sdf -i SDFID:ID -i SDFNOPREPEND -v .... file.sdf
```
`fileconv` has a very large number of options for dealing with sdf files.

All LillyMol tools support a `-i ICTE` directive, which stands for Igonore
Connection Table Errors. It will be case-by-case how to deal with molecules
that fail - ignore them, find fileconv options that attempt to fix the
problem. See [fileconv](docs/Molecule_Tools/fileconv.md).

Once you have one or more smiles files there are many things commonly done.

fileconv -v -a file.smi

provides a quick way of reporting what is in a file of molecules. For example
when run on a file derived from Chembl version 33 the output might look like
```
Molecule_Tools/fileconv_opts.cc compiled 2024-Sep-15 git hash 9654363
No output, just audit input
Processing '/home/ian/CHEMBL/chembl_33.smi'

2372674 molecules had between 1 and 780 atoms. Average 30.7147
```
Nothing was done with the molecules, smiles were interpreted and
summary statistics gathered. That took 7.5 seconds.

If you find there are very large, or very small molecules in the
collection, you may wish to discard those. 
```
fileconv -E autocreate -e -V -I 0 -O def -c 10 -C 40 -f lod -g all -v -S ok collection.smi
```
Non periodic table elements are recognised by `-E autocreate` but then discarded
by use of the `-e` option.

Molecules with likely valence errors are discarded with the `-V` option.

Molecules containing isotopic atoms are discarded with the `-I 0` option.

Molecules containing non organic atoms are discarded with the `-O def` option.
Note that molecules are first stripped to the largest fragment `-f lod`
so molecules with a Sodium counterion will be OK. 

Filters to molecules between 10 and 40 heavy atoms, `-c 10 -C 40`.

Chemical standardisation is then used to transform the structures to LillyMol
standard forms `-g all`.

Finally the filtered set is written to `ok.smi`.

Because the `-v` option is used, when done fileconv will report summary statistics
of what it just did. By default, LillyMol tools work silently.

If you have a large set of molecules, and want to impose more constraints, you
might want to use [molecule_filter](/docs/Molecule_Tools/molecule_filter.md)
instead. Filtering is an inherently parallel task, see also [dopattern](/docs/General/dopattern.md).

Once you have the set of molecules you will use for subsequent processing there
are some common tasks.

## Generate Molecular Property Profiles
One of the best ways to understand what is in a collection of molecules is to
examine profiles of common molecular properties. Use
[Molecular_Property_Profile](/contrib/Molecular_Property_Profile/README.md)
for generating property profiles.

## Uniqueness Determinations
Use either [unique_molecules](/docs/Molecule_Tools/unique_molecules.md) or
[common_names](/docs/Molecule_Tools/common_names.md) to understand
the extent of duplicates in the collection.

## MedChem Rules
run `contrib/bin/Lilly_Medchem_Rules.sh` to understand the extent to which
this collection might be fit for purpose for biological applications.

## Overlap With Existing Collections.
If you have existing collections of molecules, you can build BerkeleyDB
databases of unique smiles [structure_databases](/docs/Molecule_Tools/structure_database.md)
for those collections, and lookup the new collection in those existing
collections. You may want to discard molecules that are duplicates of
things already available.

```
in_database_bdb -d corporate.bdb -d chembl.bdb -d vendor.bdb ...
                -F duplicate -U new -g all -c -l -v file.smi
```
will look up the molecules in `file.smi` in each of those databases
and write the molecules that are NOT in them to `new.smi`.

## Diversity
There are a great many ways in which the diversity of a set can be
assessed. One simple measure is the number of unique scaffolds, from
which can get what might be called a 'scaffold density', or the average
number of molecules per scaffold.

Use [molecular_abstraction](/docs/Molecule_Tools/molecular_abstraction.md)
to convert molecules to their scaffold, and then use a text based LillyMol
tool to count the number of unique items in column 1 of the input.
```
molecular_abstraction -a 'scaffold(WRITE)' chembl.smi | unique_rows -v -c 1 -
```
Takes about 65 seconds and yields
```
Read 2236247, 1535757 duplicates, wrote 700490
```
So among 2.23M Chembl molecules there are 700490 different scaffolds
so about 3.2 examples per scaffold on average.

See fingerprinting [gfp](/docs/GFP/Diversity.md) for how to best
assess internal and inter-collection diversity.

## Substructure Searching
[tsubstructure](/docs/Molecule_Tools/tsubstructure.md) can be used
for substructure searching, and dividing molecules into those that
match a query and those which do not.

## Sorting
It can be useful and or interesting, to sort a set of molecules by
some property. Use [msort](/docs/Molecule_Tools/msort.md) for
sorting files of molecules.

## Reactions
Making changes to molecules can be done via serveral tools including
* molecular_abstraction
* [trxn](/docs/Molecule_Tools/trxn.md) 
* [minor_changes](/docs/Molecule_Tools/minor_changes.md)

## Descriptor Computation
See [/docs/Molecule_Tools/descriptors.md]

Many of these molecular descriptors are useful in QSAR model building.

## Similarity
See [gfp](/docs/GFP/AAReadme.md).

