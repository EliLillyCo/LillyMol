# Fingerprints as Descriptors

## Background

A common need is to explore use of gfp fingerprints in other contexts, for
example with similarity metrics other than Tanimoto, or to use the features
in a model that requires tabular data input.

There are two approaches to this.

## gfp_to_descriptors gfp_to_descriptors_multiple

Both of these tools read a fingerprint file and convert the fingerprint 
data contained into tabular form. The first was developed first and can
only handle a single fingerprint. The second one was built later and can
handle multiple fingerprints in the input file. Generally the second one
should be used.

### TLDR
```
gfp_make.sh -ECC3 rand.smi > rand.gfp
gfp_to_descriptors_multiple rand.gfp > rand.dat
```
For 2000 random Chembl molecules, this generates an output fle with 2001 rows
(includes a header) and 30.8k columns. This means that across those 2000
molecules there wee 30.8k different features generated.

### Usage
The usage message of `gfp_to_descriptors_multiple` is 
```
Converts a fingerprint set from a gfp file to a descriptor file
 -F <tag>       specify tag of the fingerprint to be processed, standard .gfp syntax
 -x <fraction>  discard bits that hit less than <fraction> of the time
 -n <number>    discard bits hit in less than <number> of the fingerprints
 -X <fraction>  discard bits that hit more than <fraction> of the time
 -N <number>    discard bits hit in more than <number> of the fingerprints
                    use either -x or -n, and either -X or -N
 -y <number>    discard records unless they have at least <number> non zero values
 -s             sort columns by frequency (within fingeprint)
 -I <string>    prefix for descriptors produced
 -u             gsub space to _ in multi-token names
 -d <nbits>     fold sparse fingerprints to constant width <nbits>
 -q             if a bit is absent write as -1 rather than 0
 -o <sep>       output separator (def ' ')
 -v             verbose output
```
Just like any gfp tool, you don't need to specify the fingerprint (-F)
unless you with to.

### Support
Many of the features generated will appear in only one, or a small number of
molecules, and may not be useful. You can specify the support that a feature
must have for inclusion via one of two different means.

The `-n -N` combination allows specifying lower and upper bounds on the number
of molecules that contain a feature. Just as a very rare feature may not be
informative, neither will an extremely common feature. Alterntatively the
support requirements can be expressed in fractional form via the `-x -X` options.

In the example above running
```
gfp_to_descriptors_multiple -n 10 -N 1990 rand.gfp 
```
reports
```
Will discard bits with fewer than 20 molecules hit
Will discard bits with more than 1980 molecules hit
Auto sized for
 1 sparse fingerprints
Setting 1 sparse fingerprint weights to 1
Non colliding fingerprint 0 'NCEC3C' weight 1
Read 2000 fingerprints from '/tmp/rand.gfp'
Found 30874 bits set in NCEC3C fingerprint
30177 bits below threshold 20
0 bits above threshold 1980
Will produce 697 descriptors
```
So a requirement that a feature appear in 1% of the molecules, results in
the number of features generated from 20k to 697. This dramatic reduction is
because this is a randomly chosen set of molecules. Most sets of molecules
being studied will have considerable internal similarity, and imposing
a support requirement will not have such a dramatic impact. For example
a 2750 member SAR dataset reports
```
Will discard bits set less than 0.01 of the time
Will discard bits set more than 0.99 of the time
Auto sized for
 1 sparse fingerprints
 Setting 1 sparse fingerprint weights to 1
 Non colliding fingerprint 0 'NCEC3C' weight 1
 Read 2746 fingerprints from '/tmp/data.gfp'
 Will discard bits that occur in fewer than 0.01 molecules
 Will discard bits that occur in more than 2719 molecules
 Found 15412 bits set in NCEC3C fingerprint
 14525 bits below threshold 27
 2 bits above threshold 2719
 Will produce 885 descriptors
```
so 2750 related molecules generate more features than 2000 random molecules.
And there were only 15k bits found, rather than the 30k found in the random set.

### Fixed Width Fingerprints
`gfp_to_descriptors_multiple` can also generate a fixed width tabular output.
In the case of sparse fingerprints, this will necessarily come at the cost of
possible collisions, where arbitrary feature numbers get hashed to the same
value. For example on a recent SAR dataset running
```
gfp_to_descriptors_multiple.sh -d 256 -v -x 0.01 -X 0.99 train.gfp
```
reports
```
Will discard bits set less than 0.01 of the time
Will discard bits set more than 0.99 of the time
Auto sized for
 1 sparse fingerprints
 Setting 1 sparse fingerprint weights to 1
 Non colliding fingerprint 0 'NCEC3C' weight 1
 Read 2746 fingerprints from '/tmp/data.gfp'
 Will discard bits that occur in fewer than 0.01 molecules
 Will discard bits that occur in more than 2719 molecules
 Will fold fingerprints to a constant width of 256 bits
 Found 15412 bits set in NCEC3C fingerprint
 After folding to 256 bits, how many of the initial bits hit each fixed width bit
 126 of 256 fixed width bits never set
 Of bits set, set between 1 and 483 ave 6.80769
 14525 bits below threshold 27
 2 bits above threshold 2719
 Will produce 256 descriptors
 77 of 256 bits with collisions (fraction 0.300781) max collision count 482
```
Indeed, folding those 15k bits found, or the survivors of the support requirements
to just 256 bits, does indeed result in significant bit collisions.

Note that bits are scanned before folding in order to impose support requirements.

## Generate Fixed Width Fingerprints Directly
While `gfp_to_descriptors_multiple` is a suitable means of converting fingerprints
to tabular form, some fingerprint generators can generate tabular output
directly.

### iwfp
This tool generates linear fingerprints - like the traditional Daylight fingerprint.
This is a very old tool, and internally, fingerprints are generated as fixed width -
most more recent tools use a Sparse_Fingerprint_Creator object to hold a hashed
set of bits and counts. So `iwfp` can directly generate tabular output of fingerprints.
```
iwfp -P UST:AQY -a -a -g all -l file.smi > file.dat
```
will generate a 2001*2049 file (file.smi contains 2000 molecules). The number
of columns can be controlled via the -c option - it defaults to 2048. Note that
two `-a` options are needed in order to get counted features, just one instance
results in binary output.

### iwecfp
This is LillyMol's extended connectivity (Morgan) fingerprint generator. You can
generate descriptor file output via the `-Y desc=<n>` option combination, generating
\<n\> features.

```
iwecfp -R 3 -v -P UST:AHY -Y desc=1024 file.smi > file.dat
```
generates a descriptor with 1024 columns. Note that because we used the `-v` option
we get a report on the bits set per molecule
```
Fingerprints had between 16 and 138 ave 75.996 bits set
```
which is typical of EC type fingerprints. Note that there is not a strong need
for many descriptors, 1024, of even narrower, will see few collisions.

Note too that if we use a simpler atom typing, `-P UST:Y` we get
```
Fingerprints had between 14 and 116 ave 63.9895 bits set
```
We can see the influence of using the more complex fingerprint. A very
complex atom type `-P UST:ACFHLOPSZ` results in
```
Fingerprints had between 16 and 160 ave 81.679 bits set
```
Adding one round of Morgan type expansion to this `-P UST:ACFHLOPSZ1` yields
```
Fingerprints had between 16 and 186 ave 93.869 bits set
```
Contrast this with linear fingerprints where
```
iwfp -B -P UST:ACFHLOPSZ1 -a -a
```
yields
```
Bits hit between 32 and 1741 average 756.212
```
a dramatically different result, typical of the difference in bits generated
by EC and linear type fingerprints. But generally EC type fingerprints do better
in models.


### extended_atom_pairs
This is LillyMol's atom pair fingerprint generator. You can generate descriptor
file output via the `-G desc=<n>` option combination, generating \<n\> features.

```
extended_atom_pairs -G desc=1024 -v -P UST:ACY
```
reports
```
Fingerprint 'NCYAC<' between 10 and 528 ave 179.722 bits set
```

