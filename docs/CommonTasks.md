# Common Tasks
When using LillyMol many simple tasks tend to be performed regularly. Here
is a summary of some of the most common/useful.

# Standardization
It is very important that molecules be standardized the same way. See
[standardization](/docs/Molecule_Lib/chemical_standardisation.md) for
info. Generally this means adding `-g all` to at least the first invocation
of a LillyMol tool.

## fileconv
Get an idea of what is in a file
```
fileconv -v -a file.smi
```
does no ouput, just reads the file and, because the `-v` option is used
reports a result like
```
2360 molecules had between 4 and 40 atoms. Average 20.261
```
If the file has errors, try adding `-A 2` - which allows rings
with 2 electrons to be aromatic. If there are still problems, identify
them with
```
fileconv -A D -v -a -B 999 -B log=bad file.smi
```
will ignore as many as 999 otherwise fatal errors, and write those
connection tables to the file `bad`. Examine that file to see
what the problem might be.

Note that all LillyMol tools have a `-i ICTE` option combination for
Ignoring Connection Table Errors, `fileconv` is the only tool that supports
retrieving and logging the otherwise fatal problematic inputs via the `-B` option.

Options that might be useful at this point
```
fileconv -J all -Y FHrmsqb -s rmbad -i xcte -g all ...
```
and several others. To filter to an "organic" subset, something like
```
fileconv -c 10 -C 50 -I 0 -O def -E autocreate -e -V -f lod -f rmxt=15 -v -g all -S organic ... file.smi
```
. allow creation of non periodic table
elements `-E autocreate` but then discard those molecules `-e`.

. drop molecules containing invalid valences `-V`.

. drop molecules containing isotopes `-I 0`, but maybe you want `-I change`.

. Drop molecules containing non organic atoms `-O def`. Maybe you want to consider Selenium 
organic, `-O Se`.

. If a molecule consists of two fragments each with >= 15 heavy atoms, consider it a mixture
and drop it. Otherwise reduce to the likely largest fragment.

. Apply Chemical standardisation `-g all` - this includes removing explicit Hydrogen atoms.

. Drop molecules with fewer than 10 or more than 100 heavy atoms.

. Write suviving molecules to `organic.smi`.

If you had a .sdf file as input and wanted a .sdf file as output, just add `-o sdf` and
it will write `unique.sdf` instead.

## tsubstructure
This is the first point of call for all things related to substructure searching.
The typical task is to figure out how many molecules in a given file match one
or more of a set of queries. For example if you are doing a reaction enumeration and
want to identify carboxyllic acids
```
tsubstructure -s '1[OH]-C=O' -m carboxyllic_acid -v all.smi
```
will write those molecules that have one instance of a carboxyllic acid to
`carboxyllic_acid.smi`. The `-v` option means it will report the results of
the matching.

Substructure searching can be very hard and frustrating. If you are dealing with
a complex query and having trouble figuring out what is wrong, try adding three
`-v` optins, `tsubstructure -v -v -v -s 'smarts'...` and every time it attempts
a match, it will report the number of substructure atoms that were successfully
matched. That will usually be a great pointer to where the problem lies. Note
that substructure atoms are usually matched left to right, but they depth
first follow branches.

A very common task is to leave an isotopic label on the matched atom(s). In the
case of the acid, we may want an isotope 1 on the '[OH] atom. In that case we
could either use a recursive smarts or specify that the C=O atoms are not
to be included in the resulting set of matched atoms.
```
tsubstructure -s '1[$([OH]-C=O)]' -j 1 -m carboxyllic_acid ...
tsubstructure -s '1[OH]-[/IWxC]=[/IWxO]' -j 1 -m carboxyllic_acid ...
```
each of these approaches is complex in their own way. Both yield the
same result. And of course each query could be constructed to leave the
isotope in any atom.

By default, the `-j` option to tsubstructure will increment the isotopic
number applied to matched atoms. This is seldom what is wanted, so most
often the combination `-j 1 -j same` is used. Or some other number. This
should have been the default. There is considerable flexibility in what
can be done with isotopes and queries, try `tsubstructure -j help` for
info. Each has proven useful one or more times.

## Uniqueness
This is a very common task within Cheminformatics. See the discussion
on uniqueness in [smiles](/docs/Molecule_Lib/smiles.md) on theoretical
aspects of this.

To scan through a file and just drop the duplicates
```
unique_molecules -g all -l -c -z -v -S unique file.smi
```
will be a common case, the unique molecules are written to `unique.smi`.
Without a `-S` option, `unique_molecules` just runs a diagnostic. See
[unique_molecules](/docs/Molecule_Tools/unique_molecules.md).

If you care about what is a duplicate of what
```
common_names -g all -l -c -z -v -S common file.smi
```
see [common_names](/docs/Molecule_Tools/common_names.md);

## Fragmenting
Many times there is need to decompose molecules into fragments. [dicer](/docs/Molecule_Tools/dicer.md]
is one way of doing that.

## Scaffold Hopping
Replacing the rings in a molecule is one way of doing scaffold hopping
[scaffold hopping](/docs/molecule/ring_replacement.md).

## Molecular Asbstractions
Converting a molecule into abstract forms is often an interesting way of identifying
commonality among molecules [molecular abstractions](/docs/Molecule_Tools/molecular_abstraction.md)
can be useful for that task.

## Sorting
When dealing with complex Cheminformatics tasks it is usually beneficial to start
working with smaller molecules. Sorting by heavy atom count can be done via
```
msort -k natoms file.smi > file.sorted.smi
```
will likely do what is needed. See [msort](/docs/Molecule_Tools/msort.md].

## Reactions
See [trxn](/docs/Molecule_Tools/trxn.md).
