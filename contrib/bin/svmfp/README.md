# SVMFP
This directory contains infrastructure for building SVM fingerprint models.

This requires a license for [svm_lite](https://www.cs.cornell.edu/people/tj/svm_light/)
if you are using it for commercial purposes. Contact the author of svm_lite in order
to ensure you are in compliance with the license requirements of that software.

If ok, I can supply you a customised version of svm_learn that processes a Tanimoto
kernel function and which is used by this software.

More recently, we worked with NVidia in order to make this functionality available
within their software when running on GPU's. I am however not familiar with
how this ultimately ended up in their product suite - I left Lilly just after
we got proof of concept. Presumably it is now generally available. If anyone
knows more, please let me know.

## Background
SVM is a kernel based method. Generally this means that the only thing the learner
knows about is the pairwise similarity between members of the training set.

LillyMol includes [gfp_make](/docs/GFP/gfp_make.md) which can compute an
almost unlimited number of fingerprints and fingerprint combinations.
These disparate measures of molecular similarity can be used as kernel functions.
That is the basis of SVMFP models.

Within Lilly there is a tool that explores different fingerprints, kernel functions,
across a range of train/test splits of a dataset. Work is underway to make
available similar functionality here.

## TLDR
An svmfp model requires two files, `train.smi` and `train.activity`. The smiles
file is a space separated file consisting of records like
```
smiles1 id1
smiles2 id2
.
```

while the activity file is a tabular file with a header, that must contain activity
data for every id in `train.smi`

```
ID Response
id1 1.23
id2 4.32
..
```

The ordering of the rows between the two files does not matter. Extra data in 
the activity file is ignored.

In order to build a model
```
svmfp_make.sh -mdir <dir> -A train.activity train.smi
```
will generate a default set of fingerprints, call svm_learn and write the resulting
model to `dir`.

Evaluating a model built with svmfp_make can be done via
```
svmfp_evaluate.sh -mdir <dir> test.smi > test.pred
```

## Details.
Across a variety of targets, experience tells us that the best performing fingerprint
is usually a combination of fingerprints, sometimes quite complex combinations with
3-6 components. Usually some kind of EC fingerprint will be a component, with
CATS type fingerprints often combined.

For many datasets, adding ALogp will be beneficial, especially many ADME related targets.
It is important to note that the implementation of ALogp contributes only a small number
of bits to the overall similarity - kernel function. It is often the case that the
number of bits assigned to ALogp needs to be increased beyond the default. This will
typically be done by specifying the ALogp fingerprint as `-ALOGP40` which adds
for replicates of the ALogp fingerprint - rather than 10 which is the default.

In order to build a model with a custom fingerprint combination try something like
```
svmfp_make.sh -mdir MODEL -A trian.activity -GFP -ECC3 -RS -ALOGP40 -CATS12 -GFP train.smi
```
Note both the opening and closing -GFP options. Everything between those two is
passed to gfp_make.sh. See [gfp_make](/docs/GFP/gfp_make.md) for details.


