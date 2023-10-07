Peptides as Molecules

One of the features of LillyMol is the ability to handle arbitrary strings
as elements. So within LillyMol, this is a perfectly valid smiles
```
[Th][Eq][U]IC[K][Br]O[W]NFO[Xj][U][M]PSO[Ve][Rt][He][La][Zy][D]O[G]
```
as is
```
[Ala][Arg][Asn][Asp][Cys][Gln][Glu][Gly][His][Ile][Leu][Lys][Met][Phe][Ser][Thr][Trp][Tyr][Val]
```
There are of course no concepts of aromaticity, but rings and branches are fully
supported.
```
C1C([Xy][Z])C[Qq]1
```
is a substituted 4 membered ring. Note the mixing of normal elements, C, with the
non periodic table variants. 

Note too that as new elements are discovered, the LillyMol periodic table will
be updated. For example once upon a time, 'Cn' would describe an (impossible)
atom that was 'an aliphatic Carbon and an aromatic Nitrogen'. Today that is 
element number 112, Copernicium.

These extensions are enabled via command line options

. -E autocreate

All [A-Z][a-z] two letter combinations are supported, the first example above.

. -E anylength

Atomic symbols can be any sequence of letters, so `[Serine]` would also be OK. Note
that `[serine]` is unrelated, since it has a different name.

## Unique Smiles
Identify can be done via unique smiles formation.

## Substructure Searching
Molecules like these can be substructure searched via a simple extension to smarts
```
tsubstructure -E anylength -s '[#{Ala}]!@[#{Arg}] ...
```
looks for an Alanine with a non ring bond to an Arginine.

The no matched atoms between directive can also be interesting with peptides
```
[$([ND1]-[CD3]-C(=O)[OH])]...[#{Ala}]...{>4}[#{Gly}]*[#{Ser}]
```
will find situations where starting on a terminal amine Nitrogen, any number of
atoms to an Analine, and then more than 4 atoms to a Glycine, skip one then
a Serine.

And of course ring membership can be used. This is exactly the same as normal
substructure searching, the only difference is that the elements are not limited
to two characters - and some operations involving the elements will be slower
since these are not hashed the same way that the normal periodic table is.

## Similarity Searching
Similarity searching needs to be done with a custom atom type. None of these
non periodic table elements have an atomic number, and many default atom types
involve the atomic number. As an alternative, the 'G' atom type can be used,
which uses an arbitrary hash value associated with each element - this mechanism
could be made much more efficient if there was a need.

To generate linear fingerprints, one might use
```
iwfp -R 10 -E anylength -P UST:G file.smi > file.gfp
```
which could probably be done via `gfp_make` jumping through some command line
hoops.

Atom pair fingerprints could be generated via
```
extended_atom_pairs -P UST:G -C 12 -E anylength file.smi > file.gfp
```

Once a fingerprint file exists, it is just a completely normal gfp file and
all gfp tools can be used. Experiments show that indeed very similar 'molecules'
differ by single residue changes.

And presumably too the impact of a single residue change would be more
pronounced in the middle of the 'molecule' than at the peripery.

For example looking for nearest neighbours within a file could be done with
```
gfp_nearneighbours_single_file -n 10 -v file.gfp > file.nn
```
Perhaps sort so that those 'molecules' with the closest nearest neighbours
are at the top of the list
```
tdt_sort -T DIST file.nn > file.nn.sorted
```
Convert that to smiles form
```
nplotnn file.nn.sorted > file.nn.sorted.smi
```

Or if you wanted the file in tabular form,
```
nn2csv file.nn.sorted.smi > file.nn.sorted.csv
```
which however cannot be viewed as a structure by any software.
