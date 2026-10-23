# Peptides
Over the years a common task has been to enumerate small peptides.

This tool, while quite old, can perform that task.

The input must be a file with the proposed peptide entered one sequence
per line. The names of the peptides can be in any form. Here is the
data file containing the building blocks
```
O=[1C]([1OH])C[1NH2] G Gly glycine
O=[1C]([1OH])[C@@H]([1NH2])C A Ala alanine
SC[C@H]([1NH2])[1C](=O)[1OH] C Cys cysteine
OC[C@H]([1NH2])[1C](=O)[1OH] S Ser serine
O=[1C]([1OH])[C@H]1[1NH]CCC1 P Pro proline
O[C@H](C)[C@H]([1NH2])[1C](=O)[1OH] T Thr threonine
O=[1C]([1OH])[C@@H]([1NH2])C(C)C V Val valine
O=[1C]([1OH])[C@@H]([1NH2])CCSC M Met methionine
[1NH2][C@@H](CC(=O)N)[1C](=O)[1OH] N Asn asparagine
OC(=O)C[C@H]([1NH2])[1C](=O)[1OH] D Asp aspartic_acid
O=[1C]([1OH])[C@@H]([1NH2])[C@@H](C)CC I Ile isoleucine
O=[1C]([1OH])[C@@H]([1NH2])CC(C)C L Leu leucine
OC(=O)CC[C@H]([1NH2])[1C](=O)[1OH] E Glu glutamic_acid
O=C(N)CC[C@H]([1NH2])[1C](=O)[1OH] Q Gln glutamine
O=[1C]([1OH])[C@@H]([1NH2])CCCCN K Lys lysine
O=[1C]([1OH])[C@@H]([1NH2])CC1=CN=CN1 H His histidine
O=[1C]([1OH])[C@@H]([1NH2])CCCNC(=N)N R Arg arginine
O=[1C]([1OH])[C@@H]([1NH2])CC1=CC=CC=C1 F Phe phenylalanine
OC1=CC=C(C[C@H]([1NH2])[1C](=O)[1OH])C=C1 Y Tyr tyrosine
O=[1C]([1OH])[C@@H]([1NH2])CC1=CNC2=CC=CC=C12 W Trp tryptophan
```
Names can be any of the names above, but must be consistent within
your input file.

For example to build 'P F H R', given the input file 'seq.txt'
```
P F H R
```
which would be the same as 'Pro Phe His Arg' which would be the same
as 'proline phenylalanine histidine arginine', run
```
ruby contrib/data/Peptides/peptides.rb -l seq.txt
```
The `-l` option says to use the single letter abbreviations. This generates 
```
O=[1C]([C@H]1[1NH]CCC1)N[C@H]([1C](=O)N[C@H]([1C](=O)N[C@H]([1C](=O)[1OH])CCCNC(=N)N)CC1=CN=CN1)CC1=CC=CC=C1 P-F-H-R
```
and of course any number of peptides can be specified, and any number of
lines can be in the input.

## Peptides as Atoms

As mentioned elsewhere, arbitrary symbols can be used as 'atoms' in LillyMol
so 
```
[G][A][C][S][P][T][V][M][N][D][I][L][E][Q][K][H][R][F][Y][W]
```
together with `-E autocreate` is a valid molecule in LillyMol.
Unfortunately, since some of the letters are associated with actual elements for
which there is valence information, you will get lots of warning
messges about valence errors. And since some of the letters are associated
with smarts directives you will not be able to do substructure
searches. On the other hand
```
[Gly][Ala][Asp][Cys][Ser][Ile][Leu][Glu][Gln][Lys][His][Arg][Phe][Tyr][Trp]
```
together with '-E autocreate -E anylength` will not generate any valence
errors, and can be used in substructure searches. To look for a two connected
Aspartic Acid `[#{Asp}D2]` wwill work. `[DD2]` is an invalid smarts. 

Matched with some of the range based queries in LillyMol this enables
queries like
```
[#{Ser}]...{8-12}[#{Ala}][#{Leu}][#{Glu}]
```
which is a Serine, between 8 and 12 atoms to a sequence of [Ala][Leu][Glu].

And of course normal atoms can be interspersed, the `Ala` are just atoms (Elements)
for which there is no information about atomic number, valence, mass...
```
[Ala]NC(=O)(Cc1ccccc1)[Leu][Gly]
```

is handled properly as a smiles, although whether this is useful or not is an open question.

Note too that these "molecules" can be fingerprinted. Atomic properties
referred to as "atomic number" actually use a unique hash associated with each
element. So even though 'Ala' does not have an atomic number, it does
have a unique identifier and that is what is used for an atom type. Many
other molecular properties will be applicable: connections, ring membership,
etc.
Most atomic properties will clearly not be applicable: aromaticity, Hydrogen
count, unsaturation.
