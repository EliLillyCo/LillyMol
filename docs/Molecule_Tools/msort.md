# msort

## Purpose

`msort` is a tool that sorts a file of molecules based on a range
of simple molecular properties.

## HOWTO
The usage message for msort is
```
Sorts a structure file by various criteria
Molecule_Tools/msort.cc compiled 2022-Oct-10 git hash 0e61add
  -k <key(s)>    sort specification(s)
                 a preceeding negative sign indicates reverse order for that property
                 'a' 'natoms' number of atoms
                 'r' 'nring'  number of rings
                 'w' 'amw'    molecular weight
                     'amwnH'  molecular weight, excluding Hydrogens - groups tautomers
                 'f' 'nfrag'  number of fragments
                 'c' 'nchiral' number of explicit chiral centres
                 'h' 'hetero' number of heteroatoms
                 'j' 'aroma'  number of aromatic atoms
                 'k' 'aromr'  number of aromatic rings
                 'R' 'lgrsz'  atoms in largest ring
                 'S' 'lgrss'  rings in largest ring system
                 'b' 'rotbond' number rotatable bonds
                     'asr'    atoms in smallest ring
                     'alrss'  atoms in largest ring system
                     'unsat' number of (non aromatic) unsaturated bonds
                     'ailf'   atoms in largest fragment
                     'aiff'   atoms in first fragment
                     'amwlf'  amw in largest fragment
                     'aicm=<n>' atoms in counterions, missing assigned <n>
                     'amwcm=<x>' amw in counterions, missing assigned <x>
                     'col=nn' numeric contents of column <nn> of the name
                     'niso' number of isotopic atoms
                     'rngat' number of ring atoms
                     'Z' sum of atomic numbers in the molecule
                     'sp3' number of sp3 atoms
                     'nbonds' number of bonds (single=1, double=2, triple=3)
                     'charge' number atoms with formal charges
                     'qry=...' hits to substructure query
                     'smt=...' hits to substructure query
                     'sdf=...' values in SDF tag
  -d             descending order
  -y             if multiple queries present, treat as a group
  -D <stem>      write different groups to output files starting with <stem>
  -e <number>    hint for minimum number of molecules per output file
  -M ...         miscellaneous other options, enter '-M help' for info
  -q             quick exit - avoids overhead for deallocation
  -A ...         Aromaticity, enter "-A help" for options
  -E ...         Standard element creation options.
  -i <type>      specify input file type. Enter '-i help' for details
  -v             verbose output
```

The most common atrribute is typically the number of atoms. A typical invocation 
might be
```
msort -l -k natoms -k amw file.smi > file.sorted.smi
```
This sorts `file.smi` by atom count, within the largest frgament since
the `-l` option was specified. A secondary sort criterion is molecular weight,
so this will place the Fluroine version of molecule before the Bromine.

When `msort` development commenced, sorting options were specified by
letter options. As more sort criteria were added, a more flexible way was
needed, and so the -k option was used. That should be used going forward.

It is important to note that `msort` differs from other tools in an
important way. The result, the sorted molecules, is a string echo of
the starting molecule. What gets stored is a file offset and byte 
count. When writing happens, it seeks to the appropriate offset, and
writes the byte count. Because of this it can sort .sdf files or
similar. Storing only the sort keys and file offset, also lowers
the memory footprint. But it also means that the resulting molecules
will be unchanged.

## Attributes
The following sorting attributes are implemented, and can be specified with the
`-k` option.

### -k natoms
The number of atoms in the connection table. If there are explicit Hydrogen atoms
those are also counted.

### -k nring
The number of rings in the molecule.

### -k amw
The average molecular weight of the whole molecule.

### -k amwnH
Molecular weight, excluding Hydrogens, can help group tautomers together.

### -k nfrag
The number of fragments in the molecule.

### -k nchiral
The number of explicit chiral centres in the molecule.

### -k hetero
The number of heteratoms in the molecule. Could also be specified
via a smarts, but this is more efficient. Several properties are like
this.

### -k aroma
The number of aromatic atoms in the molecule.

### -k aromr
The number of aromatic rings in the molecule.

### -k lgrsz
The size of the largest ring in the molecule.

### -k lgrss
Number of rings in the largest ring system. Benzene->1, naphthane->2.

### -k rotbond
Number of rotatable bonds.

### -k asr
Atoms in the smallest ring (size of smallest ring).

### -k alrss
Atoms in largest ring system. Benzene->6, naphthane->10.

### -k unsat
Number of (non aromatic) unsaturated bonds.

### -k ailf
atoms in largest fragment.

### -k aiff
Atoms in first fragment.

### -k amwlf
amw in largest fragment.

### -k aicm=\<n\>
Atoms in counterions, if there is no counterion, arbitrarily assign
\<n\> to the count.

### -k amwcm=\<x\>
amw in counterions, missing assigned \<x\>.

### -k col=nn
Numeric contents of column \<nn\> of the name.

### -k niso
Number of isotopic atoms.

### -k rngat
Number of ring atoms.

### -k Z
sum of atomic numbers in the molecule.

### -k sp3
number of sp3 atoms.

### -k nbonds
number of bonds (single=1, double=2, triple=3)

### -k charge
number atoms with formal charges

### -k qry=...
Number of hits to substructure query. Note that multiple
matches are counted, so if you are looking for presence or
absence of a feature, you will need to limit that in the query.

### -k smt=...
Number of hits to substructure query. Same considerations as the
`-l qry` option.

### -k sdf=...
values in SDF tag. Will need to use `-i info` in order
to have this information available.

## -d
Reverse the sort, descending order. But note that individual
sort keys can be reversed as needed.

## -y
Governs how are multiple queries handled. By default, the sort
criterion is number of matches to each query, handled separately.
So, if there are two queries, there are two sort criteria. With the
`-y` option, all queries are combined into a single sort criterion,
and the value for sorting is the sum of the hits to each one.

## -D \<stem\>
This can prove useful when dealing with a large collection of
molecules. Sorting molecules is generally a fast process - assuming
that the molecular property being perceived are cheap. If you want
to find duplicates within a large set of molecules, sorting that set
first and then looking, in parallel, at different chunks can be
a good approach.

For example, imagine a large set of molecules where duplicates are
to be identified. First sort that set via some criteria. Here is an
example involving a 21.5M set. The objective is to divide the set
into shards containing approximately 200k molecules each. The important
thing is that molecules in different shards are guaranteed to be
different from each other.
```
msort.sh -k natoms -k amw -k nring -k aroma -D shard -e 200000  -v file.smi > file.sorted.smi
```
This generates 55 files, and takes under 7 minutes, while consuming 1.3GB of RAM
on old hardware. While there may be duplicates within a shard
there are no duplicates across shards. `unique_molecules` can then
be used to process shards in parallel, although it should be noted
that the shards with small molecules are processed much quicker
than the shards containing larger molecules.

There may be other reasons for wanting disparate subsets of molecules.

## -e \<number\>
A hint on the approximate number of molecules in each shard.

## -M svsf
A variant on the `-D` splitting is to enforce having only records
that compare equal being in the same file. What is different
is that each file contains as many records as there are molecules with
that single property value, whereas with the default `-D` handling, there
might be multiple different values in each file - but again all
values that compare equal are guaranteed to be in the same file.

If the sort key is the number of atoms, the first split file contains
those molecules with the lowest atom count.

The abbreviation 'svsf' stands for Same Values Same File.

## -M firstprop
This works in conjunction with the `svsf` directive. With this active,
the file name index of the files produced is the first sort property.
So, if you want to generate a bunch of files grouped by heavy atom
count, that might look like
```
msort -k natoms -D atoms -M svsf -M firstprop file.smi > file.sorted.smi
```
which will generate as many `atom%d.smi` files as there are different
numbers of heavy atoms in the set. Note that there may be gaps
in this sequence, and it will not include zero.

## -q
Do not use. Historically, deallocating all the memory allocated
by the process was slow. The `-q` option does a quick exit 
instead of waiting. No longer an issue.

## -M
Various other behaviour modifiers are available via the `-M` option. Some
have been described above as part of the explanation for file chunking.

### -M ea=\<n\>
In order to deal with the problem of shards containing molecules of
diferent size, this option makes the shards contain approximately
equal number of **atoms**, summed across all molecules in the shard.

### -M s2f
The `col=nn` directive can be used to sort numeric tokens in the
molecule name. But if there is a text column, that cannot be
sorted that way. Instead, text is transformed into an arbitrary
number, and those numbers sorted. For example `foo` might be
translated into 1, and `bar` into 2, which would mean that all
occurrences of `foo` would occur before all occurrences of `bar`.

### -M rpt=\<n\>
Report progress every \<n\> molecules read. This should probably
be a `-r` option to main.

## `msort_parallel`
There is a parallel version of msort that processes the example
above in 2 and a half minutes. `msort_parallel` is mostly compatible
with the serial version - the code bases should be merged, but
queries are not allowed in the parallel version since currently
the substructure query object is not thread safe - it stores 
information about the current match within itself.
