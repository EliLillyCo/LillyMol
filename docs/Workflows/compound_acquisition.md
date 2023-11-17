# Compound Acquisition

Some thoughts on and methods for compound acquisiton.

## Background
For many years Lilly maintained an active compound acquisition effort, with
the objective of adding pharmaceutically reasonable, diverse molecules to the
screening set. That effort used LillyMol tools.

Note that this writeup will refer to some Lilly specific commands, for example `in_lilly_database.sh`.
In most cases it should be obvious what that command is doing. There will also be
examples of commands that relate to the specific on-prem cluster environment. Again
the meaning should be obvious.

## Mechanics
Many vendor offerings will come as '.sdf' files. The first task is to transform
that file to smiles. This can be very challenging. If you have ever wondered why
there are so many options to fileconv, dealing with external '.sdf' files is one
of the reasons. A minimal invocation might look like
```
fileconv -c 15 -C 50 -E autocreate -f lod -O def -V -e -i SDFID:CATNUM -i sdf
        -i ignore_bad_m -i mdlquiet -i ICTE -I 0 -i xctb -v supplier.sdf
```
and possibly several other options.

Adjust the min (-c) and max (-C) heavy atom counts to taste. There is no right answer.

This will create `supplier.smi`.

Discard any we already own
```
in_lilly_database -c 1 -l 50 -v -U not_owned supplier.smi
```
which creates `not_owned.smi` containing those molecules we do not already own.

At this stage it may be useful to generate molecular property profiles of the
proposed set of molecules - in order to gain understanding of their overall
characteristics. Bearing in mind that we have purposefully truncated by heavy
atom count in the previous step.

### Medchem Rules
Run Lilly Medchem Rules in order to discard molecules containing undesirable
structural motifs.

```
medchem_rules not_owned.smi > not_owned.okmedchem.smi
```
Adjust the soft and hard cutoffs in the rules to taste. Remember there is a hard
lower atom count cutoff, but soft and hard upper atom count limits.

Possibly summarise the results to see if there are any particular rules having
major impact on shrinkage. Also examine survivors to see if the rules need
to be updated for motifs not previously encountered.

Note that if this is not done, since the current set has been curated with
respect to the Lilly Medchem Rules, augmenting it with unfiltered molecules
would presumably see a strong tendency to fill those missing motifs.

### Desirability
This step is optional. Skip if you want to have a selection guided by
diversity only.

At this stage there are a range of things that could be done to assign relative
desirability scores to the candidate molecules, and or do further filtering.

* clogp, rotbond, TPSA or other models (ADME...)
* ring rarity - favor molecules with new/rare rings
* perhaps discard demerited molecules altogether, or do a demerit based scaling.
* other ideas

There will of course be no 'correct' means of combining different utility
functions, but there is extensive experience within Lilly on this subject.

### Compare vs Lilly
Generate fingerprints and split
```
gfp_make.sh not_owned.okmedchem.smi > not_owned.okmedchem.gfp
iwsplit -tdt -n 1000 -suffix gfp not_owned.okmedchem.gfp
```
which generates a bunch of iwsplit*.gfp files. Note that 1000 items in each
split results in a run-time of about 10 minutes.

Use the cluster to compare these split files with the inventory set.
```
parallel_nn_search.sh -n 1  -v -split iwsplit
```
If you need to run against some other set of molecules, the diversity
cassette, get a fingerprint file for that set, gzip it and provide it
as an argument 
```
gfp_make.sh diversity.smi > diversity.gfp
gzip diversity.gfp
parallel_nn_search.sh -n 1 -split iwsplit -haystack diversity.gfp.gz
```

this runs asynchronously and will create a bunch of files like `NNPS26337q.e128891346.5`
which are arbitrary names generated by the queuing system. 
To avoid the arbitrary names, and to get synchronous execution, you might
prefer the `-S` option.

These can be gathered into a form useful for gfp_spread via
```
parallel_nn_search_to_gfp_spread.sh -F not_owned.okmedchem.gfp NNPS26337q.o128891346.* > for_spread.gfp
```

If you look inside that file it might look like
```
$SMI<ClC1=C2C(=CC(=C1Cl)OCC(=O)O)CC(C2=O)CC>
... fingerprints omitted
PCN<PBCHM33794>
NNSMI<lilly smiles 1>
NNID<lilly id 1>
NNDIST<0.17>
|
$SMI<S(SCCNCC(C1CCCCC1)CC)CCNCC(C1CCCCC1)CC>
... fingerprints omitted
PCN<PBCHM38167>
NNSMI<lilly smiles 2>
NNID<lilly id 2>
NNDIST<0.31>
|
```
where we see that each Pubchem molecule has been annotated with the smiles, LSN
and distance of its nearest neighbour in the Lilly collection. This information,
the 'NNDIST' value will be given to gfp_spread to indicate that this is the
starting distance for this molecule. The molecule with the longest distance from
a Lilly molecule will be the first selected. A candidate molecule that is
very close to a Lilly molecule will appear non diverse and will not be
an early selection.

### Run Spread
If you have a large collection of molecules use of a parallel version of spread
is recommended. Usually spread does not benefit from more than 8 cores
```
gfp_spread_standard.sh -v -N NNDIST -v -h 8 -r 1000 for_spread.gfp > for_spread.spr
```
which may take considerable time to run depending. You may want to only generate
a number of selections that is similar to your budget - but beware that it is well
known that spread will preferentially select weird and strange (aka diverse) molecules
first so if you need 1000 molecules, you might need to ask spread to identify 20k
molecules. Doing the spread calculation is expensive - post-processing the results
in various ways is cheap.

The `-r 1000` option is just for monitoring. For every 1000 clusters formed
report progress to stderr. Gives you an idea how long this might take.

The other thing that may substantially lower the computational cost is to
decide that there is simply no way we want any molecule that is closer
than (say) 0.2 to an existing Lilly molecule. Filter those out of the
for_spread.spr file. Not sure there is a script to do this, so it may
involve multiple steps. But this step could very significantly cut down
the size of the spread run. If you have 4M candidates and can only select 1000
it probably does not make sense to do a calculation on more than the top
100k most diverse candidates.

### nplotnn
Post-process `for_spread.spr` with `nplotnn`, examining the distance associated
with each selection. Make sure that the distance to a previously selected molecule
remains desirable. It could be there may not be a large number of desirable molecules
among the candidates.

Again, the first molecules from spread may be quite undesirable - spread will
correctly identify the most diverse molecules. You may end up having to discard
the first few thousand selections depending...

### Desirability Weighted
If you can come up with a relative desirability for each molecule, best is to
have a number in (0, 1] and add that to each molecule's name. So 
```
C methane
```
turns into
```
C methane 0.7
```
which when fingeprinted will look like
```
PCN<methane 0.7>
```
With the relative desirability encoded this way, spread can be run again, but this time
adding
```
-P COL=2
```
which means that the relative desirability of each molecule is the second token
in the name field. As each distance is computed, it will be multipled by this
factor, so something with a low weight will think that it is not very diverse,
and will have a lower probability of being selected than a very similar molecule
with a higher weight. A weighting scheme like this can help 'tame' the
weird and strange molecules that spread generates, but likely you will still
need to do some triage.

A relative desirability weighted approach does seem desirable.


# Summary
The method outlined here was used for many years and still seems to be a
reasonable approach. It very directly addresses the issue of diversity
selecting from the candidates those that maximally complement what
is already here. The cluster is leveraged for the inevitably large
pairwise comparison of each candidate to the existing collection.
Multi-threaded processing is used for large spread runs - the maximum
diversity selection, although in most circumstances that job could
be trimmed to moderate sizes with low risk. This would lead to
very attractive run-times.

Using external desirability measures enables further refinement of
selections towards a mixture of desirability and diversity.

The main disadvantage of the method is the tendency of spread to
(correctly) identify weird and strange molecules as diverse, and
generally some level of postprocessing is needed in order to suppress
these.