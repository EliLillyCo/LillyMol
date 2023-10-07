# A population for training a deep learning model.

Many deep learning generative models need to be trained on a
large corpus of valid molecules. Preferably molecules that have
actually been made, since this should maximise the likelihood
that the model will then generate plausible molecules.

Because of the desire to process exemplified molecules, collections
like the Lilly collection and Chembl are to be preferred over
things like Pubchem, Enamine Real or Zinc, each of which has
a great many virtual molecules. It is quite possible that
many of those molecules could be made, but there are too
many of them anyway...

## Do Not Want - Simple Filters
For our purposes, we will restrict attention to molecules containing only
the 'organic' subset of elements, 'C, N, O, F, P, S, Cl, Br, I'. We choose to 
exclude elements such as B, Si and Se that are sometimes considered.

We restrict attention to molecules having between 10 and 50 heavy atoms, since
it is very likely that all atomic arrangments of interest are exemplified
in molecules of that size range. If we were to change this range, perhaps
[15,40] would be a better choice.

We exclude molecules with apparent valence errors, isotopes and only consider
the largest fragment.

## Halogens
It is extremely rare for a QSAR model to be improved by differentiating the
heavy halogens, so we commonly consider all heavy halogens to be equivalent.
Many LillyMol tools have an element transformation option, (typically `-t` or `-T`)
which specify element transformations. That is done here. Usually we
consider Cl more desirable than Br -> more desirable than I, so the usual
invocation might be
```
-t I=Cl -t Br=Cl
```
which converts all I and Br to Cl.

## Smiles
If the learner is learning from smiles strings, it would seem likely that
the presence of Chlorine might make the learning process harder - since the
model will need to learn that 'Cl' is different from 'C'. For that reason,
we replace all heavy halogens with 'I' since it is a single letter. We are
not looking to generate molecules containing Iodine, but when molecules are
generated 'I' atoms can be replaced with 'Cl'.
```
-t Cl=I -t Br=I
```
will do that. Note that changing an existing unique smiles via text manipulation
may not work.

We may choose to discard molecules containing too many halogens. For example
from a 100k random sample from Pubchem we find this many 'Iodine' atoms
```
 29101 molecules had 1 hits
 5106 molecules had 2 hits
 1283 molecules had 3 hits
 247 molecules had 4 hits
 42 molecules had 5 hits
 21 molecules had 6 hits
 9 molecules had 7 hits
 3 molecules had 8 hits
 1 molecules had 9 hits
 1 molecules had 10 hits
 1 molecules had 12 hits
 1 molecules had 13 hits
 1 molecules had 14 hits
```
This includes things like
```
ClC12C3C4C(C5(Cl)C(Cl)(Cl)C4(Cl)C(Cl)=C5Cl)C4=CC(=C(S(=O)(=O)O)C=C4C3C(Cl)(C1(Cl)Cl)C(Cl)=C2Cl)C(=O)O PBCHM5771665
```
which has 12 actual Chlorine atoms!

In order to speed processing, it may be convenient to use 'grep' to remove
molecules with too many 'Iodine' atoms, or `tsubstructure -s '>3I' -n ok -m toomanyhalogen ...`.
We really don't want our generative model thinking that molecules with 12 heavy halogen
atoms are OK. Remember, we have not considered F atoms here, just the heavy halogens.

### Aromatic Nitrogens
In an aromatic smiles, aromatic Nitrogen atoms are marked with the number
of implicit Hydrogens - otherwise the smiles can be ambiguous across tautomeric forms.
Smiles will contain a mixture of '[n]' and '[nH]' symbols. Perhaps it might be
desirable to remove this complexity, as long as there are tools to read
those possibly ambiguous smiles. This remains to be tested.

## Unique Smiles
It is unclear whether training a model on canonical smiles is beneficial or not. We
should examine the efficacy of learning a represention on

- 1M different, unique smiles
- 500k different molecules, two non unique smiles of each
- Other??

## Chirality
Many parmaceutically relevant molecules contain chirality markers.

| Collection | Fraction |
| ---- | ---- |
| LLY | 0.21 |
| Chembl | 0.24 |
| Pubchem | 0.22 |

a surprisingly consistent fraction across disparate collections. However
the quality of this information is very uncertain. Some of us take a quite
skeptical view of chirality information and would generally prefer to discard
it. While we have a pretty good idea of the accuracy of the chirality information
in the corporate collection, we assume that the external collections are
even less precise.

There may be a case to be made for retaining chirality if a generative model
is to be built with 3D targets in mind, although even then, just trying all
plausible enatiomeric forms would seem worth trying.

### Filters
Applying molecular desirability filters can be expensive. Applying the
Lilly Medchem Rules to 100k random Pubchem molecules takes 15 seconds, or about
2.5 minutes per million. If many millions of molecules are to be processed
parallel processing should be used. 

I was going to recommend that the the `-nodemerit` option be used, but
that would suppress things like the Nitro rule, which allows one, but not
more instances of a Nitro group. We do not want our learner to think that
molecules with 3 Nitro groups are OK.

Filters based on ML models may be too slow to be practical, although
even if a subset of the molecules were processed, that might be
beneficial.

# HowTo

If you are going to filter with the Lilly Medchem Rules, that may need to
be done first. The Halogen transformations described above will violate
certain rules, depending on how many Halogens you retained.
You can also impose size limits here if required.

The following protocol seems to work - medchem rules done later...

For large collections from which a subset is needed, use `random_records`
to extract a random set. We anticipate this being faster than a combination 
of 'shuf' and 'tail', but not sure. Certainly less memory intensive.

Once the subset is identified, filter with fileconv, enforcing the various filters
and transformations described previously.
```
fileconv -c 15 -C 40 -O def -t Cl=I -t Br=I -K nonH -I 0 -f lod -s 0 -V -E autocreate -g all -v -o usmi -S collection.norm collection.smi
```
should do that, generating 'collection.norm.smi' from 'collection.smi'. Running this
on 100k random Pubchem molecules takes about 5 seconds, so less than a minute
per million.

If this is done for each collection, `unique_rows` can be used to eliminate
the duplicates.

```
unique_rows -c 1 -v collection1.norm.smi collection2.norm.smi ... > unique.smi
```
should do that, or to handle excess heavy halogens
```
grep -v 'I.*I.*I' collection*norm.smi | unique_rows -c 1 -v - > unique.smi
```

# Summary
Using LillyMol tools to assemble a set of molecules useful for training a
large language model seems feasible.

