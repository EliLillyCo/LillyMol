# Synthetic Precedent

LillyMol contains a tool for doing synthetic feasibility assessments. This was
independently developed at about the same time that Peter Ertl implemented a
similar tool in RDKit.

The LillyMol tool does not do any special processing of rings or ring systems,
the [ring rarity](smi2rings.md) tool is used for handling rings.

The first step is to build a database of EC fingerprints from molecules
which you believe have actually been made. From these molecules we 
generate EC fingerprints, and store the info in a BerkeleyDB database.
Importantly, we also store a count of the number of instances of that
bit found in the collections. In many larger collections, we find
molecules that are likely to have been drawing errors, or otherwise
mysterious looking motifs. A simple presence or absence of a feature is generally
not a robust determinant of rarity.

## Building A Database
The tool iwecfp_database_load builds a database from one or more
collections. These should be molecules where there is reasonable
expectations that the molecules have actually been made.
The keys in the database are bit numbers and
the values contain information about prevalence, and maybe examplar structures.
Newer versions will contain protos.

A typical database load might look like
```
iwecfp_database_load -d precedent.bdb -R 3 -P SFX collection.smi
```
For a collection like Chembl, this might take 20 minutes.

## Lookups
Once the database has been built, new molecules assessed.

A typical lookup command might look like
```
iwecfp_database_lookup -d precedent.bdb -w PSD unknown.smi
```
If we build the database with Chembl, and then lookup Chembl in that
database, the output might look like
```
N1=C(C2=CC=CC=C2N=C1N)NCC CHEMBL4218505 45 bits score 0.1 0 302759 1 8123 2 432 3 1
O1C2=CC(=CC=C2OC1)OC1=CC=C(N)C=C1 CHEMBL1382750 50 bits score 1.6232 0 180412 1 9792 2 1106 3 42
ClC1=C(S(=O)(=O)N2CCCCC2)C=CC=C1Cl CHEMBL1525354 46 bits score 1.415 0 263575 1 14656 2 79 3 26
O=C(C1=CC=C(O)C=C1O)C(=O)CCCCC CHEMBL1579361 54 bits score 0.1 0 640344 1 3018 2 1 3 1
S1C(=NCC1)NC(=O)C12CC3CC(C1)CC(C2)C3 CHEMBL1310254 46 bits score 0.1 0 29815 1 1474 2 4 3 1
C1=[N+](C=CC2=C1C(=CC=C2)OC)C1=CC(=CC=C1)OC CHEMBL4784748 63 bits score 0.30103 0 9687 1 137 2 22 3 2
FC(F)CN1N=CC(=C1)N[C@@H](C)C1=NC2=C(O1)C=CC=C2 CHEMBL4914797 70 bits score 0.30103 0 90661 1 966 2 6 3 2
C1[C@@]2(C)NC3=NC4=C(N3[C@@H]1C1=CC=CC=C1O2)C=CC=C4 CHEMBL1171309 67 bits score 0.1 0 44492 1 281 2 1 3 1
C12(C(=O)N(C)C(=O)C3=C1C=C(C(F)(F)F)C=C3)C(=O)NC(=O)C2 CHEMBL294313 74 bits score 0.30103 0 100348 1 445 2 16 3 2
FC1=CC=CC=C1CCNC1=C(C(=C(N(=O)=O)C(=N1)C)C)C#N CHEMBL1381951 73 bits score 0.60206 0 109719 1 5268 2 23 3 4
```
The first column is the smiles and then the ID of the molecule being assessed.

Then follows the number of EC bits found in this molecule. 

Following that is a score. Various attempts at creating single numeric scores for this
tool have all proven unsatisfactory in one way or another. The score shown is just the
logarithm of the radius 3 count. Using the count itself seems more straightforward.
Newer versions will deprecate with this.

Then follows three sets of 'radius count' pairs.

For example for the last molecule in the table above CHEMBL1381951 the pairs are
```
0 109719 1 5268 2 23 3 4
```
This means that at radius 0, the rarest radius 0 shell (which is an atom) has
109719 examples in the database. At radius 1, the rarest shell has 5268 examples.
Radius 2 there are 23 examples, and finally at radius 3 there are 4 molecules
in the database that have this radius 3 arangement of atoms. 

A radius of 3 is quite large, seven atoms across, so this is quite specific. It is
quite possible for a molecule to have a missing bit at radius 3 and be a quite
reasonable molecule, just that particular arrangement of atoms was not in
the reference collection.

We can ask the question which are the most precedented molecules - those
with the highest number of examples at radius 3.
```
CCCCCCCCCCCCCCCC CHEMBL134994 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCCCCCCCCC CHEMBL1234557 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCCCCCCCC CHEMBL135488 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCCCCCCC CHEMBL135694 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCCCCCC CHEMBL30959 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCCCCC CHEMBL132474 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCCCC CHEMBL134537 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCCC CHEMBL335900 14 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
CCCCCCCC CHEMBL134886 13 bits score 4.419 0 1522691 1 132180 2 45790 3 26243
```
which is hardly surprising. Every molecule in the collection that has a C7 chain
will set that bit

And which molecules have the smallest number of radius 1 examples
```
S1(=O)[C@@H]([C@H](C)[C@@H](C)[C@@H]1C=CC)[SH+]O CHEMBL1173558 43 bits score 0.1 0 1 1 1 2 1 3 1
O=[N+]1C(N)=C(NCO)C(=N[C-]1N)N CHEMBL1741970 47 bits score 0.1 0 1 1 1 2 1 3 1
N1([C@H]2N=CN[C+]2[C-](N)N=C1)CC=CC1=CC=CC=C1 CHEMBL1741808 62 bits score 0.1 0 1 1 1 2 1 3 1
P1=C([N+]2=C(C1)C=CC(=C2)CCCC)C(=O)OCC CHEMBL1992180 64 bits score 0.1 0 1 1 1 2 1 3 1
P1(C2=C(C3=C1C=CC=C3)C=CC=C2)CC1(COC1)CP(CC1=CC=CC(=C1)C)CC1=CC(=CC=C1)C CHEMBL1993304 70 bits score 0.1 0 1 1 1 2 1 3 1
C1=CC(=CC(=C1OCC(O)CNC(C)COP)C)NC(=O)NCC CHEMBL3247824 81 bits score 0.1 0 1 1 1 2 1 3 1
ClC(Cl)(Cl)C1N(C2=CC=CC=C2)C2=S(N1C(=O)C)N(C1=CC=CC=C1)C(=N2)C1=CC=C(OC)C=C1 CHEMBL4889750 87 bits score 0.1 0 1 1 1 2 1 3 1
COC1=CN=C(C=N1)C(=O)NC1=CC=C(F)C(=C1)[C@]1(C)CS(=O)(=O)C2(C[O+](C)C2)C(=N)N1 CHEMBL4110849 102 bits score 0.1 0 1 1 1 2 1 3 1
P(=C(F)N(C)C)C(F)(F)F CHEMBL1980015 27 bits score 0.1 0 2 1 1 2 1 3 1
ClC(Cl)=PC1=C(C(C)(C)C)C=C(C=C1C(C)(C)C)C(C)(C)C CHEMBL1999422 36 bits score 0.1 0 2 1 1 2 1 3 1
[Mg]12(OC3=CC=CC=C3C(=O)O1)OC1=CC=CC=C1C(=O)O2 CHEMBL3561635 36 bits score 0.30103 0 2 1 2 2 2 3 2
[Mg]12(OC3=CC=CC=C3C(=O)O1)OC1=CC=CC=C1C(=O)O2 CHEMBL3580437 36 bits score 0.30103 0 2 1 2 2 2 3 2
I1(=O)(O)OC(=O)C2=CC=CC=C12 CHEMBL118857 40 bits score 0.1 0 2 1 2 2 1 3 1
```
These are some rather horrible looking molecules, possibly some may be erroneously
drawn, hard to know. For the first one, there is an atom in that molecule where all
shells have just one example in the database - this molecule itself.

We see the importance of using counts, rather than presence or absence, of EC fingerprints.
If just presence or absence is known, unknown molecules containing these features would
be assessed as OK.

## Dicer
The tool described above is the most commonly used synthetic precedent assessment tool
in LillyMol, but it is also possible to use dicer fragments.

Again, this is a two phase process where one or more existing collections
are diced, a database is built, and then new molecules can also be diced
and the resulting fragments looked up in the database.

This has the advantage of using fragments that are likely human interpretable, 
rather than EC shells, but suffers from the disadvantage of being more expensive to
compute.

### Dicing the Collection and Building the Database.
Run dicer on the collection(s). For example that might look like
```
dicer -I 1 -c -k 3 -X 500 -M 12 -B nbamide -B brcb -B fragstatbinproto \
        -B fragstat=chembl.fragstat.tfdata -T I=Cl -T Br=Cl -B nooutput -v chembl.smi
```
This may take a while. Note that the maximum fragment size, -M, is an important
consideration. Only frgaments within this range will be stored in the database. If
the molecules to be examined contain large fragments that cannot be decomposed by
dicer, those fragments will not influence the outcome. Make the -M option as
large as you reasonably can, but larger values will make the datbase (much) larger
and lookups will also be slower as more fragments need to be checked.

Other options can be adjusted as needed (dicer)(dicer.md).

The output is a TFDataRecord file of serialized protos of the fragment
decompositions found during dicing. 

That file can then be loaded into a BerkeleyDB database
```
dicer2bdb -v -d chembl.dc.bdb chembl.fragstat.tfdata
```
which forms 'chembl.dc.bdb` which can be used for lookups.

### Dicing New Molecules and Lookups.
Once the database is build, run dicer on new molecules, using the same
options that were used when building the database.
```
dicer -I 1 -c -k 3 -X 500 -M 12 -B nbamide -B brcb -B proto -T I=Cl -T Br=Cl -v new.smi > new.dc.textproto
```
this generates textproto output which can then be consumed by the database lookup tool
```
dicer_fragment_lookup_bdb -v -d chembl.dc.bdb  -t new.dc.textproto
```
output might look like
```
Cc1ccc2[nH]ccc2c1 CHEMBL112462 1328 [nH]1c2c(cc1)c[1cH]cc2 9
Oc1o[n]c2c1cccc2 CHEMBL282723 19 o1[n]c2c([1cH]1)cccc2 9
CN1CCCC1=NC1CC1 CHEMBL3303624 2 [1NH]1C(=NC2CC2)CCC1 9
Nc1[nH]c(CCNC)[n][n]1 CHEMBL1624638 2 CNCCc1[nH][1cH][n][n]1 9
CCC(C)N1CCNCC1 CHEMBL1621578 6349 C[1CH2]CC 4
O=C(O)C1CNCCC1O CHEMBL16226 25 OC(=O)C1[1CH2]CCNC1 9
NCCCC(N)(C)C(=O)O CHEMBL3247534 7 C[1CH](N)CCCN 7
C#CC(N)CC(CN)(F)F CHEMBL1178865 2 [1CH3]C(C[1CH2]C#C)(F)F 8
N#Cc1sc(N(=O)=O)cc1 CHEMBL2323851 92 N#Cc1s[1cH]cc1 7
NS(=O)(=O)OC1CCCC1 CHEMBL150300 1701 [1OH]S(=O)(=O)N 5
NS(=O)(=O)OCC(F)(F)F CHEMBL153068 755 NS(=O)(=O)O[1CH3] 6
O=c1oc2cccc(C)c2[nH]1 CHEMBL320850 59 O=c1oc2c([nH]1)[1cH]ccc2 10
O=C(C)Oc1ccc(C)cc1 CHEMBL501246 82 O=[1CH]Oc1ccc(C)cc1 10
```
The first token is the smiles of the input molecule, then the id.
Then follows information about the rarest fragment found from
that molecule. In the case of the last molecule above, CHEMBL501246,
there were 82 occurrences of fragment `O=[1CH]Oc1ccc(C)cc1` found
in the database. The trailing 10 is the number of atoms in the fragment.

The dicing and lookups can of course be pipelined
```
dicer -I 1 -c -k 3 -X 500 -M 12 -B nbamide -B brcb -B proto -T I=Cl -T Br=Cl -v new.smi | \
        dicer_fragment_lookup_bdb -v -d chembl.dc.bdb  -t -
```

At one extreme of output we have molecules that frgament into very
common fragments
```
O=C1C(CCC2OC2CCCC1(C)C)(C)C CHEMBL2000792 1701930 [1CH4] 1
O=c1[n](c2c([n]3c1[n][n](c3=O)C)cccc2)C CHEMBL421427 1701930 [1CH4] 1
O=C1NCC2CCCCNc3[n]c4c(c5[nH]c2c1c5)cccc4[n]c3C CHEMBL3808669 1701930 [1CH4] 1
O=C1N(CCC2(CC3N(CCc4c3oc3ccccc34)CC2)N1C)C CHEMBL1394278 1701930 [1CH4] 1
O=C1NC(=CCCC(=CCC(C=C1)(C)C)C)C CHEMBL3628870 1701930 [1CH4] 1
O=C1N(CC(N1C)C12CC3CC(CC(C3)C1)C2)C CHEMBL398055 1701930 [1CH4] 1
O=N1=Cc2c(CC1(C)C)c1c(cc2)cccc1 CHEMBL148184 1701930 [1CH4] 1
```
and at the other end of the range we have molecules where the rarest fragment
is in fact very rare in the database.
```
c1cc2c(cc1c1[n]c3[n]c[nH][n]3c1)CCCC2 CHEMBL566731 1 [nH]1[n]2c([n][1cH]c2)[n]c1 8
c1ccc2c(c1)CN(CCN1Cc3ccccc3OC1)CO2 CHEMBL1098598 1 [1CH3]CN1COc2c(C1)cccc2 12
CC1(c2cc(NS(=O)(=O)C)ccc2)C2CN(CC3Cc4ccccc4C3)CC12 CHEMBL1957717 1 [1CH3]N1CC2[1CH](C)C2C1 8
CC1(C)CC2(CCC(C1)(OO2)c1ccccc1)c1ccccc1 CHEMBL103236 1 C[1CH]1C[1CH]2OO[1CH](C1)CC2 10
Cc1cc(c2[n]c(OCCOC)cc(COC)[n]2)ccc1 CHEMBL1890788 1 [1CH3]COc1cc([1CH3])[n][1cH][n]1 10
CC1CCc2sc(NC(=O)c3sccc3)c(c12)C(=O)N(CC)CC CHEMBL374609 1 CC1c2c(s[1cH][1cH]2)CC1 9
Cc1ccc(C2=NNC3(CCOC(C3)(C)C)S2)cc1 CHEMBL1398613 1 S1C2(NN=[1CH]1)C[1CH2]OCC2 10
Cc1cc(ccc1S(=O)(=O)c1ccc(N)cc1)N CHEMBL174208 1 Nc1cc(C)c([1SH](=O)=O)cc1 11
```
In these cases, the rarest fragment has just one instance in the database - and
given that this database was built from Chembl, it must have been this starting
molecule.

Long term, it remains to be seen whether the EC type precedent tool, or the dicer
based precedent tool proves most useful.

