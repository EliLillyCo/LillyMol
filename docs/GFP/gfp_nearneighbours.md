# gfp_nearneighbours & gfp_lnearneighbours

These tools are the primary nearest neighbour tools in LillyMol. They
both perform the task of comparing every fingerprint in one file, the
needles, against every fingerprint in another file, the haystack. The
result is
```
for every needle
  a sorted list of neighbours from the haystack
```

Output is a TDT file that can be processed by any of the tdt_* utilities
or [nplotnn](nplotnn.md).

## TLDR
For each needle, find the 10 nearest neighbours in the corporate
collection using the default fingerprints. Sort the results so
the fingerprints with the clostest nearest neighbour are at the 
top of the list.
```
gfp_make.sh needles.smi > needles.gfp
gfp_make.sh collection.smi > collection.gfp
gfp_lnearneirhbours -p needles.gfp -n 10 collection.gfp > needles.nn
tdt_sort -T DIST needles.nn > needles.nn.sorted
nplotnn -v needles.nn.sorted > needles.nn.sorted.smi
```

Nearest neighbour finding is an intrinsically parallelisable task.
A very common approach is to split the needles, search each chunk
of needles against a larger collection and then collate the results.
```
iwsplit -n 1000 -suffix smi needles.smi
dopattern.sh -parallel 10 -suffix smi 'gfp_make.sh %.smi > %.gfp'
dopattern.sh -suffix smi -cluster 'gfp_lnearneirhbours -p %.gfp -n 1000 haystack.gfp.gz > %.nn
nplotnn -v *.nn > needles.nn
```
The smiles files are split into 1000 unit chunks. Local resources
are used (-parallel 10) to generate fingerprints for each chunk.
Then each chunk is sent to the cluster for the nearest neighbour
determination.

Note that we show a gzip'd haystack file. Whether it is beneficial to
compress the haystack or not will depend on the local computing environment.
Experiment and see.

## Two Tools?
When originally envisaged the thought was that the haystack would be
slurped into RAM, and once that was done, the needles could, one by one,
rapidly be compared with every member of the haystack. This is `gfp_nearneighbours`.

Somewhat unexpectedly, we found that it was faster to do things the
other way round. Read the comparatively small number of needles into RAM.
Then the haystack is read one fingerprint at a time and compared with each
member of the needles set. Presumably this was due to the overhead associated
with a large number of memory allocations required to slurp potentially
large numbers of fingerprints into RAM. It also imposes limitations due
to RAM availability - which were a lot more severe when these tools were
developed. This was of performing the search is `gfp_lnearneirhbours`
which tends to be more commonly used.

For that reason, the discussion will focus on that tool. Generally
`gfp_nearneighbours` is very similar wrt options, inputs and outputs.

## Options.
The most important options are
### -n \<nbrs\>
For each member of the needles, the number of neighbours to accumulate.
Run times go up as larger and larger sorted lists of neighbours need
to be maintained.

### -T \<dist\>
An upper distance cutoff. If a distance from needle to haystack member
is found to be greater than \<dist\> the fingerprint is discarded.
Note that this may mean that even if you have requested \<nbrs\>
neighbours be found, some fingerprints may have fewer than that
number of neighbours - even zero.

But there is also a -m option, which forces a minimum number
of neighbours to be found when the -T option is specified.
```
gfp_lnearneirhbours -T 0.30 -m 5 -p needles.gfp haystack.gfp > needles.nn
```
generally discards haystack fingerprints greater than 0.30 away from the
needle, but it will retain the 5 closest ones found - which will
by definition be further away than 0.30.

This raises an important point of strategy with nearest neighbour
searches. It is expensive to perform the search, it is cheap to
post-process the results with `nplotnn`. If you wish to impose
an upper distance threshold, you migh be best searching for (say)
1000 nearest neighbours (-n 1000) and then post-processing the result
with nplotnn's -T option.
```
gfp_lnearneirhbours -p needles.gfp -n 1000 haystack.gfp > needles.nn
nplotnn -v -T 0.25 needles.nn > needles.smi
```
this is generally very fast since it is just doing text processing.

Finding nearest neighbours can often be a complex undertaking.

### -t \<dist\>
This imposes a lower distance cutoff. For example you may not be
interested in duplicates so `-t 0.001` should accomplish that.
Note there is also a `-h` option. This discards neighbours that
have both zero distance __and__ the same identifier. But note
that `nlotnn` also has a -t option, so close neighbours can
also be handled during postprocessing.

### -g \<dist\>
This is an optimisation. Just as most GFP tools have a -W option,
which says do __not__ compare molecules if the atoms counts differ
by a given amount, the -g option does something similar during
multi-fingerprint distance computations. A distance computation
generally involves evaluating multiple fingerprints, and reporting
the, possibly weighted, average of those fingerprints. The -g option
says that if any of the individual distance computations is further
than the cutoff, abandon the computation at that point.

There is risk involved with this. You might miss a neighbour
if the first fingerprint evaluated showed a large distance,
but the next fingerprint would have shown a much shorter
distance. There are no free lunches here.

### -r \<n\>
Finding nearest neighbours can be a lengthy process, and
nothing is written to stdout during processing. If you are curious about
how long a task might take add `-r nnn` and every `nnn` fingerprints
processed it will report to stderr.

### -k
While the traditional TDT output generally works well within the
LillyMol infrastructure, processing via protos offers some 
advantages. That output might look like
```
smiles: "C(Cl)(Cl)(Cl)CCl" name: "CHEMBL155816"
        nbr { smi: "C(Cl)(Cl)(Cl)C(Cl)(Cl)Cl" id: "CHEMBL160929" dist: 0.13132441 }
```
which is fairly straightforward for a single nearest neighbour, but becomes
more complex when multiple neighbours are processed. More tools are
being developed to consume these, and other proto type files.

## Other options.
The options above account for almost all modern day use of the tool. 
Other options are a mixture of obsolete, bad ideas or very
specialised uses.

Many such options are available via the -B option. Invoking with `-B help`
yields
```
The following -B qualifiers are recognised.
 -B nofatal         ignore otherwise fatal errors.
 -B nosmiles        discard neighbour smiles - things run faster and consume less memory.
 -B nbrtag=<tag>    write the number of neighbours for each target in <tag>.
 -B bignn           use algorithm optimum for large neighour lists.
 -B rmzero          remove leading zero's from identifiers.
 -B ewt             distance metric is equal weight Tanimoto.
 -B maxdistanceinclusive    distances up to and including the max distance are considered matches.
 -B ID=<tag>        in the output, the tag that gets echo'd as the identifier (default PCN).
 -B DIST=<tag>      in the output, write distances with the <tag> (default DIST).
```

### -B nofatal
Occasionally distances outside the range [0,1] may result. As the code
base has become more robust, this happens very rarely today. Generally
if an out of range distance is encoutered the programme terminates immediately.
With this setting, such distances are truncated.

### -B nosmiles
Writing the smiles of each neighbour in the output file can lead to
such output files becoming quite large. With this setting only the name
of the neighbour, and the distance, are written.

### -B nbrtag=\<tag\>
By default the neighbours are written with the default identifier tag, PCN. That
can be changed with this setting. This is never used these days.

### -B bignn
If large numbers of neighbours for each needle are requested - hundreds or more,
this changes how the accumulating list of neighbours is stored - sorting at
the end, rather than attempting to maintain a sorted list throughout. The
only way to determine if this is beneficial is to try it with and without.
Over time develop a feeling for when it is useful with the kinds of searches
being done. The difference is usually small.

### -B rmzero
Strip leading zeros from identifiers. Do not use this.

### -B ewt
The distance function used is an Equal Weight Tanimoto of all bits from
all fingerprints. By default, each fingerprint is evaluated separately
and the results combined. With this setting there is a single Tanimoto
computation done, involving all the bits from each of the fingerprints.
This mimics what happens in svm_learn, where bits lose knowledge of
from where they came.

### -B maxdistanceinclusive
Change the test for maximum distance from `< maxdist` to `<= maxdist`.

### -B ID=\<tag\> -B DIST=\<tag\>
Do not use. Controls the tags used to write the output. The defaults are
always OK these days.

## Other obscure and obsolete command line options.

### -s \<size\>
This is obsolete and is generally not used.

When processing commences, the needles file is scanned to see how
many fingerprints are in there. The -s option pre-specifies the
size. You can specify a number smaller than what is in the file
and then only that many fingerprints are read.

It is also possible to have a per molecule maximum neighbour
distance via `-T TAG=<tag>` but this is seldom used. Input might
look like
```
$SMI<smiles>
PCN<id>
MAXT<0.2>
FP...
```
and invocation would be
```
gfp_lnearneirhbours -T TAG=MAXT ...
```
It is unclear that this is useful.

# Summary
gfp_lnearneirhbours is used extensively within LillyMol enviroments
for finding nearest neighbours. It is fast and adapts readily to parallel
and distributed envornments.
