# grep_molecule
This tool is inspired by the simple question
```
Is this molecule (needle) in this file of molecules (haystack).
```
for arbitrary needles and haystack files.

## Expensive
Obviously one way to do this would be to use buildsmidb_bdb to build
a BerkeleyDB datbase of the unique smiles of the 'haystack' and then
use in_database_bdb to lookup the needles in the haystack database and
report the matches.

But this is a complex, multi-step process, and far out of whack with the
simplicity of the question - a pair of files, or a smiles string and a
smiles file. Note that grep does not work unless smiles are unique.

## grep_molecule
One obvious way of solving this problem would be to
1. Read the haystack molecules and store the unique smiles in a hash.
2. Read the needles, compute their unique smiles and look in the hash.

While this would be fine for small sets of molecules, it is very
undesirable for larger haystacks. This is buildsmidb_bdb except
that the hash is never written to disk.

By making some simplifying assumptions, we can make this process more
efficient than what the naieve implementation does.

###
The most obvious optimisation is to record, for each needle, the
atom count. As haystack molecules are read, unless there is a needle
with that many atoms, just ignore the haystack molecule.

It is important to note that the atom count can be determined as
a text operation on the smiles - without ever parsing the smiles
as a Molecule. In additon to determining the atom count, the
function can also assemble a molecular formula - without Hydrogen.
These are all simple string operations and are fast.

The next optimisation is around the molecular formula. For each needle
compute the molecular formula and place that in a hash. As haystack
molecules are read, their molecular formula is also determined. If the
molecular formula is not in the molecular formula hash, ignore the
haystack molecule.

Only once a haystack molecule has passed all these pre-filtering steps
is the smiles interpreted as a Molecule. At this point, a more
precise molecular formula can be computed, which includes
Hydrogen counts, and aromaticity. This generates molecular
formula like `C8H2c6H3O2H` where aromatic and aliphatic atoms
are differentiated. This uses the `formula_distinguishing_aromatic`
Molecule method. This is also loaded into a hash.

When a member of the haystack is checked, if the earlier filters have not
rejected it, this formula variant is checked, and if there is
not a match, the molecule is rejected.

Only after all these earlier steps is a unique smiles generated, and
looked up in the unique smiles hash.

But clearly, if you are querying a set of molecules for which there already
exists a BerkeleyDB database built with buildsmidb_bdb, use that always!

### Results
Asking whether or not methane is in Chembl can be performed by
```
grep_molecule 'C methane' chembl.smi
```
takes 1.3 seconds and yields
```
C CHEMBL17564 methane
C CHEMBL2106049 methane
```
the name of the input molecule is appended to whatever was in the
haystack file.

Multiple needle molecules can be concatenated.
```
grep_molecule 'C methane,CC ethane,CCC propane,C1CC1 cyclopropane' chembl.smi
```
again takes 1.3 seconds and yields
```
C CHEMBL17564 methane
CC CHEMBL135626 ethane
C(C)C CHEMBL135416 propane
C1CC1 CHEMBL1796999 cyclopropane
C CHEMBL2106049 methane
```
The order of the output is the same as the order in the haystack file.

If this is run with the -v, verbose, option the invocation reports
```
0 failed formula generation
2409092 no atom count match
168 no formula match
4 no aromatic formula match
6 unique smiles comparisons
5 matches found
```
The haystack had 2.4M molecules. All smiles were interprted by the
string based formula generator.

Almost all the molecules, 2.4M failed the atom count filter - in other words
there are hardly any molecules in Chembl with 1-3 atoms. Of course!

Clearly this is a very artificial, and best case scenario for the omptimisations
deployed.

For something more realistic, if we take 2000 random Chembl molecules and
look for those in Chembl, that takes 20 seconds, reporting
```
Read 2286887 molecules
0 failed formula generation
783 no atom count match
1574045 no formula match
702704 no aromatic formula match
9355 unique smiles comparisons
2474 matches found
```
We see that 1.57M failed the formula match. A further 702k failed the
aromatic formula match, leaving just 9.3k unique smiles comparisons.

Worst case would be finding a set of molecules in itself, which would require
the full unique smiles computation. Nevertheless, doing this on 20k random
molecules takes just 2.7 seconds.

## Further Optimisation
The number of chiral centres could be included with the aromatic molecular formula,
or perhaps discerned from the starting smiles. If chirality is being considered
this might help. But across Chembl only 550k molecules have one or more chiral
centre, so the benefit may be limited.

If it were known that aromatised smiles were being read, the initial molecular
formula could be made more precise, but generally we disfavour storing
unique smiles.

It is possible that the number of rings could be discarned from the starting
smiles, and that could provide a benefit.

Perhaps the aromatic molecular formula could be made more precise by 
also differentiating ring aliphatic atoms.

