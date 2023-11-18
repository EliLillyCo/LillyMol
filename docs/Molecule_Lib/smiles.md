# Smiles

Smiles in LillyMol are mostly standard, and can be consumed by many
other software packages.

## Aromaticity
Different software tools may make differing assessments of what is
aromatic. When presented with aromatic smiles from package `A`, package
`B` may be unable to parse the smiles, because package `B` do not believe
that this ring can be aromatic.

Also, fundamental to LillyMol is the need for a Kekule structure, so when
aromatic smiles are read, considerable effort must be undertaken in order
to find a valid Kekule form - which may occasionally fail.

We therefore prefer to store smiles as Kekule forms. That way we maximise
the probability that any other package can read the smiles. The resulting
extra size of the smiles is more than compensated for by the faster
processing enabled.

## Standardisation
There is a separate file on standardisation in LillyMol (link when done).
One of the decisions made in LillyMol is to prefer five valent Nitrogen atoms
over charge separated representations. Neither is correct, and they each
have advantages and disadvantages. Other packages may not be able to
consume five valent Nitrogen atoms. `fileconv -g rvnv5` will attempt
to convert all five valent Nitrogen atoms into charge separated forms
for use by third party software.

## Geometry
Coordinates can be added to a LillyMol smiles. For example benzene might be
```
C{{-0.0167,1.3781,0.0096}}1=C{{0.0021,-0.0041,0.002}}C{{1.2084,-0.6789,-0.0117}}=C{{2.3961,0.0285,-0.0201}}C{{2.3773,1.4107,-0.013}}=C{{1.1709,2.0855,0.0021}}1 benzene
```
where the atomic coordiantes are appended after each atom. This is usually
smaller than a .sdf file, and retains the convenience of having each structure
on one line. Tools that just read the connection table consume these smiles
just as if the coordinates were not there.

Over the years I have devoted considerable effort to try and find a compact
representation of the coodinates, with the objective of being able to write
the coordinates as something like `C{{121823}}=C{{121904}}` but such
efforts have generally not been successful, resulting in files that
are not that much smaller than the smiles file above. There are many
problems with this idea, and it has been abandoned, even though aspects
remain in the code. I have also tried working with binary proto files,
but again, the advantages are minor.

## Chemaxon Extensions
Chemaxon provides a wide variety of extensions to smiles via the `|...|` construct
after the smiles. LillyMol will recognise some of these, but most relate
to concepts that are not implemented in LillyMol - advanced stereochemistry, or
things that LillyMol can do via other means - arbitrary atomic symbols. Generally
support for Chemaxon extensions is slim... It will recognise coordinates.

## Identity
Determining molecular identity is one of the most common tasks for
someone working with molecules. One of the first decisions that you
need to determine is at what level of granularity do you need to
determine aromaticity? For example, do 'counterions' matter? Does
chirality matter? What about tautomeric representations? Each
of these is very complex, with no obviously right answers. It all
depends....

Generally in drug discovery, a unique smiles representation is adequate
for identity determination. It is also important to acknowledge that
default smiles is deficient as a molecular representation, and many
important aspects of molecular identity cannot be encoded in a default
smiles. Chemaxon has tried to address this with their smiles extensions.
InChI provides quite robust identity matching, but at the expense of
computational cost.

For LillyMol we deliberately adopt a simplified approach, generally
assuming that the unique smiles is an adequate measure of equivalence.

The LillyMol unique smiles suffers from cases where the same
molecule may generate different unique smiles. Generally these are
unusual structures where there can be pathological interactions between
the determination of the smallest set of smallest rings and the
determination of aromaticity. Fixing these cases would be hard and
would likely slow down all calculations.

In the case of the "organic" subset of Chembl, 5-50 heavy atoms, there
is just one failure

!(CHEMBL1998572)[Images/tsmiles_failure.png]

which is largely a result of aromaticity determination in a complex
fused system. Outside the organic subset
```
119164 molecules had between 1 and 780 atoms. Average 82.5219
```
there are 52 failing molecules.
```
52 molecules had between 65 and 136 atoms. Average 80.9038
```
So across 119k molecules the failure rate is 0.04%, but when
counted against all of Chemb, fraction failing is very low, and
restricted to molecules seldom of interest 

To replicate this
```
tsmiles -v -p 10 -a -t 0 -h -w 5 -m 5 -j -u -y -q -R -f -r 10000 -F failed -A 2 chembl.smi
```
Here are some of the failures
```
C123C4(C5=C6C7=C1C1C8C9=C7C7=C%10C%11=C%12C(=C67)C6=C5C5=C7C%13=C6C%12=C6C%12=C%13C%13=C%14C%15=C%12C%12C%16C%17%18C%19=C%20C(=C8C8=C%19C%16=C(C%10=C89)C%11=C6%12)C6=C1C2=C1C2=C6C%20C(C%15%17CN(CCOCCOCCN)C%18)C%14=C2C(=C7%13)C1=C45)CN(C3)CCOCCOCCN CHEMBL1185065 aromatic smiles mismatch
C123C4(C5=C6C7=C1C1C8C9=C7C7=C%10C%11=C%12C(=C67)C6=C5C5=C7C%13=C6C%12=C6C%12=C%13C%13=C%14C%15=C%12C%12C%16C%17%18C%19=C%20C(=C8C8=C%19C%16=C(C%10=C89)C%11=C6%12)C6=C1C2=C1C2=C6C%20C(C%15%17CN(CCOCCOCCN(CC(=O)O)CC(=O)O)C%18)C%14=C2C(=C7%13)C1=C45)CN(C3)CCOCCOCCN(CC(=O)O)CC(=O)O CHEMBL407797 aromatic smiles mismatch
C123C4(C5=C6C7=C1C1=C8C9=C%10C%11=C%12C%13=C%14C(=C6C%12=C79)C6=C5C5=C7C9=C6C%14=C6C%12C9C9=C%14C%15=C%16C%17=C%18C(=C1C2=C%17C(=C45)C%16=C79)C1=C8C%10=C2C4=C%11C%13=C6C5=C4C4=C(C%15=C%18C1=C24)C%14=C5%12)CC(OC(=O)CCC(=O)O)CC3 CHEMBL216867 aromatic smiles mismatch
```
all 52 seem to be BuckyBall type structures, again seldom of interest to
Drug Discovery.

*Note* canonicalization of cis/trans bonds in LillyMol does not
work. Work started, but once we came to the realization that most of
the cis trans bonding information we had was likely incorrect, the
effort was abandoned. Simple cases mostly work, but more complex
cases do not.

This leads to an important discussion about chirality in smiles.
Ideally all molecules would be properly annotated with chirality
unambiguously marked. The reality is different. In Drug Discovery
many molecules are never resolved. Some registration systems may
force an arbitrary choice of chirality in cases of racemic
or unknown chirality. There are stereo concepts that are difficult
to deal with in smiles - a 70/30 mixture, if site 1 is R then
site 2 is S. These are difficult to deal with in computational chemistry
in general, and likely need to be enumerated.

Many collections contain atoms that are marked chiral, but which are not,
and atoms which are indeed chiral, but which have no marked atoms. For
example in the "organic" subset of Chembl, 5-50 heavy atoms we find
```
unmarked 2253510 values between 0 and 19 ave 0.291271
```
So there is an average of 0.29 unmarked chiral centres per molecule. In
tabular form
| Number Unmarked | Frequency |
| --------------- | --------- |
| 0 | 1771711 |
| 1 | 373456 |
| 2 | 75180 |
| 3 | 16767 |
| 4 | 9226 |
| 5 | 3610 |
| 6 | 1608 |
| 7 | 581 |
| 8 | 497 |
| 9 | 203 |
| 10 | 305 |
| 11 | 120 |
| 12 | 56 |
| 13 | 53 |
| 14 | 42 |
| 15 | 70 |
| 16 | 12 |
| 17 | 6 |
| 18 | 3 |
| 19 | 4 |

Similar results are found with other collections. What this means is that
any attempt to deal with chirality, starting with something like Chembl,
has to deal with a significant number of "missing values". In addition, one
must ask, of the molecules marked as being chiral, how many of those 
annotations are actually correct? Generally we take a skeptical view of
chirality and frequently it is best handled by enumerating the
different forms, see `id_chirality`.
