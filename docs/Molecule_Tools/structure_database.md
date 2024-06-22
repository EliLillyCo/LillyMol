# Structure Database
A Structure Database is a key/value store where the key is a unique smiles and the value
is an identifier. These databases are used for answering the question 'is this molecule
included in a collection?' For example, get the Chembl ID for this unknown molecule.

## Background
This problem turns out to be surprisingly difficult. In the most simple form, the
unique smiles for the entire molecule is stored, counterions, chirality and when a
new molecule is looked up, both the fragment information and the chirality must match.
That is unambiguous and straightforward.

But in reality, we usually need more nuanced approaches to the structure equality
problem. For example, do we care about the counterion? Do we have enough confidence
in the chirality information to insist on there being the exact same chirality? Or
are we specifically looking for the other molecules that might have the same
chirality to compare their biological activity - see also 
[activity_consistency](/docs/Molecule_Tools/activity_consistency.md)

Chemical standardisation is included as a default part of this system, so
many tautomeric forms will be matched - unless chemical standardisation has
been turned off. But even then there may be a need to look for other likely
tautomeric matches.

Tools buildsmidb_dbd and in_database_bdb evolved over the years to address this
problem. They worked by building separate databases for
1. The exact unique smiles of all fragments, including chirality
2. The unique smiles of the largest fragment, including chirality
3. The unique smiles of the largest fragment, without chirality
4. The unique smiles variant that enables tautomer matching.

While this generally worked well, it was always deficient in various ways - 
keeping track of what was stored in which database for example. Also, since
many molecules do not have chirality, many times there was significant duplication
between #2 and #3.

## Structure Database
This system was designed to avoid some of the shortcomings of the previous system,
while being built in a way that facilitates a more modular approach, suitable for
Python bindings.

Similar to the earlier version, this uses BerkeleyDB databases as the key/value
store, although any such system would work.

By default, databases are automatically built as a pair. The first BerkeleyDB database
contains a combination of #'s 1, 2 and 3. If the database is named 'foo', two files,
'foo.smi.dbd' and 'foo.graph.bdb' will be generated. 'foo.graph.bdb' contains
#4 unique smiles. When the database 'foo' is opened
for lookup, both of those database files will be opened.

### Building
Building can be done via
```
structure_database_load -d collection /path/to/collection.smi
```
There is only one important option, the `-s` option. Many collections of molecules
contain marked chiral centres that are not actually chiral. This creates problems
when forming unique smiles. It is highly recommended that all databases be built
with the `-s` option, which removes these non-chiral chiral centres.
Again, this creates
`collection.smi.bdb` and `collection.graph.bdb`. Building a recent version
of Chembl with 2.37M molecules takes about 9 minutes.

The database contents can be listed with `iwbdb_list` or queried with `iwbdb_fetch`. 
For the `*.smi.bdb` file that might look like (key value).
```
C CHEMBL17564:CHEMBL2106049
C#C CHEMBL116336
C#CC CHEMBL116902
```
We see that there are two Chembl entries for methane. Their ids are concatenated as
the value for the `C` key.

But since we are also storing molecules that have been stripped to the largest fragment,
and molecules from which chirality has been removed, we also get records that look like
```
C(=C1CCC2=NCCCN12)c1ccccc1 F%CHEMBL3228328:CHEMBL3303375
```
in this case, CHEMBL3303375 has the unique smiles of the key, whereas CHEMBL3228328
generates the same unique smiles after stripping to the largest fragment - the 'F%'
prefix.

Molecules that become identical once chirality is discarded will also be grouped
together
```
CC(CC#CC(O)(C)C)C1=CCC2C(=CC=C3C(=C)C(CC(O)C3)O)CCCC12C @%CHEMBL35115:@%CHEMBL284121:@%CHEMBL3350835
```
when chiraliy was considered, these were different, but once stereochemistry is
dropped, they generate the same unique smiles.

And of course we can get fragment stripping and stereochemistry equivalents grouped
together
```
s1ccc(CN2CC3N(CCOC3C2)c2[n]ccc[n]2)c1 @%F%CHEMBL3503916:@%CHEMBL3555652
```
In this case the same structure results if CHEMBL3555652 loses chirality,
and CHEMBL3503916 loses both a counterion and chirality. Note that CHEMBL3503916
will also be in the database with the unique smiles generated when it
is reduced to the largest fragment - actually it will have three entries in
the database: the smiles of the original molecule, smiles with fragment(s) stripped,
smiles with chirality also lost.

The `*.graph.bdb` database is more unusual. In order to form the key, the
molecule is converted to graph form [mol2graph](/docs/Molecule_Tools/mol2graph.md).
Then a variant of the original molecular formula is appended. This is an aromatic
distinguishing molecular formula, where Toluene would have a formula `c6H5CH3`.
In addition, the SSSR ring sizes are appended to the formula. This ensures that
any matches based on this are likely to be plausible pairs.

A typical database entry might look like
```
CC1CCC(OCCC2CCC(CN3CCN(C4CCCCC4)CC3)CC2)CC1:c18H13C8N2OH20a666 CHEMBL4455261:CHEMBL4595079
```
We see that there are 18 aromatic carbon atoms, and 8 aliphatic carbon atoms ...
In this particular case, the pair are related by a counterion relationship, the
largest fragment is the same. But in other cases more interesting relationships
can emerge.

Consider this pair. 
```
BrC1=C(C)C(=C(C)C=C1C)C1=C(O)C(O)=C(C)C(=O)C1=O CHEMBL2144683
O=C1C(=C(C2=C(C)C=C(C)C(Br)=C2C)C(=O)C(=C1C)O)O CHEMBL4296991
```
![CHEMBL2144683](/docs/Molecule_Tools/Images/CHEMBL2144683.png)
![CHEMBL4296991](/docs/Molecule_Tools/Images/CHEMBL4296991.png)
Perhaps some form of chemical standardisation might be able to
bring these into concordance, or perhaps they are actually
different molecules.

The pair 
```
C1(=C)C(C(C)(C)CC(=C1)C)CC1=C(O)C=C(Br)C(=C1)OC CHEMBL48746
C1(C(=CC(=C)CC1(C)C)C)CC1=C(O)C=C(Br)C(=C1)OC CHEMBL299828
```
![CHEMBL48746](/docs/Molecule_Tools/Images/CHEMBL48746.png)
![CHEMBL299828](/docs/Molecule_Tools/Images/CHEMBL299828.png)

Represent an interesting structural grouping that would likely
be of interest in an SAR followup.


This pair
```
C(=C(O)C1=CC(Br)=CC=C1)C(=O)C(=O)O CHEMBL22197
C(=O)(C1=CC=CC(=C1)Br)C=C(O)C(=O)O CHEMBL564494
```
![CHEMBL22197](/docs/Molecule_Tools/Images/CHEMBL22197.png)
![CHEMBL564494](/docs/Molecule_Tools/Images/CHEMBL564494.png)

is an interesting case in tautomer forms. Perhaps there
should be a standardisation to being these together, but it
would likely be difficult - and probably rare.

The pair
```
C1(=C(O)OC)C(C(=C(C)N=C1C)C(=O)OC)C1=CC(Br)=CC=C1 CHEMBL307527
C1=CC(=CC(=C1)Br)C1C(=C(C)NC(=C1C(=O)OC)C)C(=O)OC CHEMBL2063819
```
is probably just a drawing error

![CHEMBL307527](/docs/Molecule_Tools/Images/CHEMBL307527.png)
![CHEMBL2063819](/docs/Molecule_Tools/Images/CHEMBL2063819.png)


The pair
```
C12=C(NN(C1=O)C1=CC=CC=C1)C1=CC(Br)=CC=C1N=C2 CHEMBL3144697
C1(=CC=C2C(=C1)C1=NN(C(=O)C1=CN2)C1=CC=CC=C1)Br CHEMBL4243764
```
is a possibly interesting pairing from an SAR perspective. Should
these be standardised to the same form?
![CHEMBL3144697](/docs/Molecule_Tools/Images/CHEMBL3144697.png)
![CHEMBL4243764](/docs/Molecule_Tools/Images/CHEMBL4243764.png)

These "tautomer" equivalence groupings can also bring together
certain groups that would never be thought of as tautomeric
variables, but might nevertheless be interesting
```
BrC1=C(C)N2C(=C(N=C2C=C1C)C1CCCCC1)CC1=CC=CC=C1 CHEMBL2142055
BrC1=C(C)N2C(=C(N=C2C=C1C)C1=CC=CC=C1)CC1CCCCC1 CHEMBL2134378
```
![CHEMBL2142055](/docs/Molecule_Tools/Images/CHEMBL2142055.png)
![CHEMBL2134378](/docs/Molecule_Tools/Images/CHEMBL2134378.png)

Which of this pair, or both, is correct?
```
C1=CC=CC2=C1C(=NO)C(=C1C(=O)NC3=C1C=CC(=C3)Br)N2 CHEMBL243555
N1C(=O)C(=C2C(=C3C=CC=CC3=N2)NO)C2=C1C=C(Br)C=C2 CHEMBL2178734
```
![CHEMBL243555](/docs/Molecule_Tools/Images/CHEMBL243555.png)
![CHEMBL2178734](/docs/Molecule_Tools/Images/CHEMBL2178734.png)

Across Chembl there are about 10k relationships like this that are
not captured by other structure groupings. Many would also be
captured by similarity searches, but the objective here is
quick matching, and potentially turning up unexpected relationships.

This might be used as the default for looking up in a corporate
database where it might be very important to capture all plausible
structure matches with an external molecule - and it would be better
to capture more matches than fewer. Again, similarity searching can
also play a useful role in this situation, but at much greater expense.

## Lookups
The purpose of a Structure Database is performing lookups. Now
that the mechanics of how the database is constructed is understood,
we can examine how lookups are done.

The tool `structure_database_lookup` looks up molecules from the command
line and there is a [python](/docs/python/structure_database.md) binding.

The usage is
```
Looks up molecules in a structure database built by structure-database_load.
 -d <dbname>    database name - generated by structure_database_load. 
                Both <dbname>.smi.bdb and <dbname>.graph.bdb database files must be present.
 -l             strip input molecule to the largest fragment before doing the lookup.
 -c             remove chirality from the input molecules before doing database lookup.
 -t             transform the input molecule to graph form before doing the database lookup.
 -F <file>      stream for molecules found in the database
 -p             append database identifier to name
 -e <sep>       search all databases, insert <sep> string between database finds, def , def |
 -U <file>      stream for molecules NOT found in the database
 -M <file>      stream for molecules with database lookup status
 -W <fname>     write per atom count lookup rates to <fname>
 -v             verbose output
```
The -d option will be the same as the -d option given to `structure_database_load`. There
can be any number of databases specified via -d options.

The strategy of the command line tool is to apply whatever transformations you specify,
`-l`, `-c`, `-t` and then lookup the resulting molecule in the database. If you specify
no options, then the standardised unique smiles is the lookup key. And if the input
molecule has neither fragments nor chirality, the options will have no effect.

The molecules that are found in the database can be written to the `-F` file. If you
also specify the `-p` option, the values from the database are appended to the
input molecule name in the `-F` file. Molecules not found inthe database can be written
to the `-U` file.

Alternatively, you can get summary lookup status for each input molecule with the `-M`
option. It is also common to find that smaller molecules are more often found than
larger molecules, so the `-W` option provides the ability to get a per-atom-count
summary of database match rate.


## Problems
Tools like this expose many of the weaknesses of specific unique smiles systems, and
this is no exception. For example

![CHEMBL313657](/docs/Molecule_Tools/Images/CHEMBL313657.png)

Is the chiral centre opposite the Nitrogen chiral or not? One of the carbon
atoms with OH attached is marked chiral, one is not, so we don't know what it is.

This one is probably clearer, since the atom that makes this asymmetric
is resolved,

![CHEMBL1590135](/docs/Molecule_Tools/Images/CHEMBL1590135.png)

The Structure Database tools make an attempt to identify possibly invalid
chirality, if the `-s` option is used, but it is imperfect. It would be
much more expensive to treat more completely.

But this also raises the question of how reliable do you believe the
chirality information in the molecules to be? In most public collections
there are significant numbers of unmarked chiral atoms - presumably reflecting
the fact that they have never been resolved, and/or they might be racemic.

Some collections may have adopted registration policies where all asymmetric
atoms must be marked, even if not known. 

With the combination of conceptual and computational difficulties, the
large amount of missing information, and the likely low reliability of
much of the specified chirality, we would generally recommend performing
matches without chirality wherever possible. Or possibly do matches with
and without and then decide whether the differences are meaningful.

If chirality information *is* considered reliable, then the unique smiles
used in LillyMol should properly handle almost all cases.

## Performance
Performance for lookup is usually similar to load performance. If the database
file(s) are "cold", things will be slow. But unlike actual database systems,
the more simultaneous users, the better the performance - since the database
files tend to get cached in RAM.

## Summary
The Structure Database tools provide a convenient and performant approach
to the structure matching problem, available from both the command line
and via python bindings.
