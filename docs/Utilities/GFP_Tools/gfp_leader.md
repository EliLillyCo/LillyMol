# gfp_leader

A clustering tool also known as sphere exclusion.

## Use Case
The most common use case for `gfp_leader` is when you have a set of molecules
which can be ordered by some desirability function, and when you need to select
a limited number of desirable, yet diverse candidates.

The desirability measure might be things like

. A docking score
. One or model model scores
. Medchem demerits
. Heavy atom count
. Other...

There is almost always some relative desirability score that can be assigned to
a set of molecules, even if just heavy atom count in order to favor smaller molecules.

## Algorithm
A radius must be specified. The guarantee is that no two selected molecules will
be closer than this distance. What distance to use will be a function of what
fingerprints are being used. For the standard `gfp` set, typical distances might
be in the range of 0.1 to 0.2, although there is no 'correct' value. It depends on
how many molecules you have, how many you must select, and how confident you are in
the scoring function.

If the scoring function is very good, you might be comfortable with a larger radius,
confident that the model is giving accurate scores to non leader molecules in
each cluster. If on the other hand you had a very weak desirability, heavy atom
count, you might want a shorter radius, since the 'model' is not very predictive
of the other molecules in each cluster. There are many factors to consider, and
choice of threshold cannot be strictly prescribed.

Some have criticized leader for the need to have an arbitrary distance. More
expensive clustering algorithms may not properly emphasise the most desirabile
molecules, and/or have their own tuning parameters.

It is also possible that any given choice of the threshold will yield too many,
or too few clusters. You may need to run multiple times with different thresholds
in order to get a desirable outcome.

Leader scales well. Worst case would be if every molecule ended up in its own
cluster, a very small radius and no duplicates. That would take `n * (n - 1) / 2`
calculations. The other extreme would be if every molecule ended up in the same
cluster, which would take `n - 1` calculations. So it is always better than `n^2`.

Leader starts with the order set of candidates. The first item is marked as selected
and then all subsequent items that are within `threshold` of that first item are
placed in the first cluster. It then returns to the beginning of the list and finds
the first unselected item. Subsequent items that are within `threshold` of that 
leader are placed in its cluster. This continues until all items have been selected,
or until the number of items needed have been selected.

## HowTo
starting with a smiles file, a likely sequence of commands might be
```
gfp_make file.smi > file.gfp
gfp_leader -t 0.15 -v file.gfp > file.ldr
```

If you have a large file, you may want to run `gfp_leader_tbb` in order to 
run multi-threaded. If you are using the standard fingerprints `gfp_make -STD` then
you can use `gfp_leader_std` which is also multi-threaded.

The following options are supported by the base, serial, version.
```
Performs leader clustering on a set of fingerprints
Usage <options> <input_file>
 -C <number>      maximum number of clusters to find
 -t <dis>         specify distance threshold
 -t col=<nn>      threshold is column <nn> of the name field
 -t tag=<TAG>     threshold for each molecule in dataitem <TAG>
 -H <TAG>         threshold for each molecule in dataitem <TAG>
 -m <number>      maximum cluster size
 -M <tag>         max cluster size for each molecule in <TAG>
 -M col=nn        max cluster size is column <nn> of the name field
 -S <TAG>         score tag
 -S col=nn        score is column <nn> in the name field
 -I <TAG>         specify identifier tag
 -r               sort clusters by distance from leader
 -E <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)
 -E ALL           echo all dataitems from the pool file
 -A <file>        file(s) of previously selected molecules - discard all within threshold
 -a <dist>        use <dist> as the threshold when comparing against the -A file
 -L <file>        write fingerprints discarded by -A file(s)
 -D ...           miscellaneous options, enter '-D help' for info
 -s <number>      specify max pool size (not needed)
 -F ...           gfp options, enter '-F help' for details
 -V ...           Tversky specification, enter '-V help' for details
 -v               verbose output
```

### -C ncluster
The number of clusters to form. By default, clustering will continue until
all fingerprints have been placed into a cluster.

### -t distance
The threshold used for cluster formation. Note that a constant threshold is
used for all selections. This may, or may not be what you want.

### -t tag=TAG
Specify a per-molecule threshold. For the most desirable molecules you may want
to use a small threshold, so only the very close neighbours get clustered together
with the most desirable molecules. Then use a larger threshold with less desirable
molecules.
```
$SMI<C>
PCN<methane>
THRESHOLD<0.11>
|
```
This is the same as `-H THRESHOLD`. It is more common to place the threshold
in a column of the name, see below.

### -t col=col
Instead of a TDT tag, per-molecule thresholds can be specified as part of the name
field
```
$SMI<C>
PCN<methane 0.11>
|
```
In the example above, `-t col=2` means that the threshold for each molecule if
found in column 2 of the name.

### -m max
Maximum cluster size. Once a cluster has \<max\> members, stop adding
items to that cluster. 

### -M TAG
Specify a per molecule maximum cluster size via a TDT tag.
```
$SMI<C>
PCN<methane>
MAX_CLUSTER_SIZE<50>
...
```
In this case, 'methane', if it is selected as a leader, can only form a
cluster with as many as 50 items. More commonly specifying a maximum
cluster size is done via a column in the name.
```
$SMI<C>
PCN<methane 50>
...
```
and then specifying `-M col=2`.

### -S tag -S col=nn
The relative score for each molecule is in a TDT tag. Usually it is easier
to sort the input file by score ahead of time.

### -A fname
Specify a file of fingerprints containing molecules that are to be treated as
having been previously selected. For example if running two different selection
methods, a docking score and a QSAR model, you can either

* Combine the scores into a composite measure and sort by that.
* Run leader with one method, then again with the other method.

A workflow that needs to select 1000 moleules might look something like.
```
sort_by_docking_score file.smi > file.sorted1.smi
gfp_make file.sorted1.smi > file.sorted1.gfp
gfp_leader -C 500 -t 0.2 -v file.sorted1.gfp > file.sorted1.ldr
nplotnn -n 0 file.sorted1.ldr > file.sorted1.sel

# The selections from the first method are now available. Fetch those
# fingerprints from the previously generated file, or regenerate.
fetch_tdt_quick -c 2 file.sorted1.sel file.sorted1.gfp > file.sorted1.sel.gfp
sort_by_qsar_score file.smi > file.sorted2.smi
gfp_make file.sorted2.smi > file.sorted2.gfp
gfp_leader -A file.sorted1.sel.gfp -C 500 -t 0.2 -v file.sorted2.ldr
nplotnn -n 0 file.sorted2.ldr file.sorted2.ldr.smi
```
The two files, 'file.sorted1.sel' and 'file.sorted2.sel' contain the
selections.

### -a threshold
By default, unselected items are discarded if they are within the
global threshold (-t) of any previously selected item (-A). That
can be changed via the `-a` option.

### -r
By default, when a cluster is written, it will be in the order of the
molecules in the input file - which was ordered by desirability. In some
cases it may be desirable to have cluster members sorted by distance
from the leader.

### -s \<size\>
The number of fingerprints in the input. Seldom needed. By default the tool counts
the fingerprints, and then allocates internal arrays. Once upon a time this
may have been slow.
