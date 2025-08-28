# Parallel Processing

LillyMol was developed in an environment with an abundance of parallel
opportunities. Many common Cheminformatics tasks are "pleasingly parallel"
because molecules can be processed independently. In those cases,
ready parallel processing is available by splitting the input info
chunks, processing each chunk independently, and then, optionally, rejoining
the result files.

A common workflow might be
```
iwsplit -n 100000 -suffix smi large.smi
dopattern -cluster -suffix smi 'expensive_computation %.smi > %.out'
```
which submits each file to a Grid Engine cluster.

Over time, the number of cores on machines has increased subtstantially, with
32 or more cores common today. For many tasks that might be adequate
```
iwsplit -n 100000 -suffix smi large.smi
dopattern -suffix smi -parallel 32 'expensive_computation %.smi > %.out 2> %.log'
```
In this case all the .smi files in the directory are processed in parallel,
with a maximum of 32 tasks active at any time.

# Parallel Tools
As noted in the explanation of IO in LillyMol [IO](/docs/Molecule_Lib/io.md)
most LillyMol Molecule aware tools support `-i seek=` and `-i stop=` directives
which enable reading from positions in a file other than the beginning.

Several tools leverage this in order to provide convenient parallel processing of files.

It is important to note that in all cases, these tools just automate the splitting
and local processing processes, although it would be easy to extend that to batch
submission.

Currently the following utilities are provided

. fileconv_parallel
. iwdescr_parallel
. tsubstructure_parallel
. trxn_parallel
. gfp_make_parallel
. unique_molecules_parallel
. get_substituents_parallel
. buildsmidb_parallel
. in_database_parallel

Almost all of them behave exactly the same as their single processor
variants - because they are calling the same executables, with options
that might look like
```
fileconv -o usmi -S tmpfile1 -i seek=1830542 -i stop= 329843319 large.smi
```
This has the very desirable feature that the input file does not need
to be split. On the other hand, the separate output files to then need
to be joined, so there is no free lunch.

Scaling is essentially linear, since the parallelism is at the process level
and the operating system is very good at managing separate processes. Given
32 cores, we have observed parallel substructure searches of Chembl taking
1-3 seconds, using the exact same options as regular tsubstructure.

gfp_make_parallel is slightly unusual because it is a perl script that
just assembles a pipelined command of different fingerprint generation
executables. It figures out which is the first component of that pipeline
and applies the seek options to just that invocation.

get_substituents_parallel is different in that each invocation generates
a dicer_fragments textproto file and these are then aggregated via
`dicer_fragments_collate`. It processes a recent Chembl in 23 seconds.

buildsmidb_parallel is also slightly different and must be used
with caution. Each task builds its own database and those databases
are then concatenated with `iwbdb_cat -a`. This is only OK that is how
you would handle duplicate entries. For example if you were using 
the -p option, do NOT store duplicates, this would not work since
duplicates could result as databases are concatenated with the -a option.
For most use cases this works fine, but be aware of these potential
pitfalls.

Note too that unlike the underlying tool, buildsmidb_bdb, this
wrapper will silently remove the target database before invoking iwbdb_cat.
This may be a bug or a feature depending on your perspective.

Takes 1 minute 17 seconds to build a BerkeleyDB database of 2.3M
Chembl structures when running with 16 threads. But note that
increasing that to 32 threads took 1 minute 14 seconds, presumably
due to I/O saturation. Beware resource limitations that may not
be immediately obvious.

in_database_parallel can be used for multiprocess database lookups.
Note that it only supports the -F and -U options for output. Looking up
1.84M Chembl molecules in the full Chembl database, takes about 6
seconds when done 16 way parallel. This time running 32 way parallel
drops the time to 3.7 seconds from which we might conclude that
I/O saturation is much less a problem with reads than writes.

## unique_molecules_parallel
This is handled quite differently. In this case the tool first uses
msort_parallel to divide the input into separate files. The sorting
is based on the molecular formula, so all molecules that have the
same formula will be in the same file. Once this is done, the task
is pleasingly parallel - unique_molecules can be invoked on each file.

Note that msort_parallel is a C++ multi-threaded application, not a
wrapper script like the tools being discussed here.

The call to msort_parallel might look like
```
msort_parallel.sh -p -D /a/temp/file -M nc=16 -k formula large.smi
```
which divides the input into 16 chunks, grouped by molecular formula.
unique_molecules is then invoked on each of the files
```
/a/temp/file1.smi
/a/temp/file2.smi
```
On a test set with 1.8M molecules, running 16 way parallel takes 19 seconds.
Running 32 way parallel takes 13 seconds - the bottleneck is msort_parallel
which is a multi-threaded tool and does not scale linearly as threads
are added.

On a test with 37M molecules, msort sorts that into 16 different chunks in
1 minute 42 seconds.

In all cases, these convenience tools are designed to lessen the time
we spend waiting for results

## Other Multi Threaded Tools
Several LillMol tools are multi-threaded. This is mostly done using OMP,
but several tools use POSIX threads.

All of these tools exhibit diminishing benefit as threads are added.

. gfp_leader_standard
. gfp_leader_tbb
. gfp_spread_standard
. gfp_nearneighbours_single_file_tbb
. gfp_lnearneighbours_standard
. nn_spread

For example, running gfp_spread_standard on a recent Chembl subset with 
about 2.34M molecules, results in
| Threads's | Time (sec) | Hours |
| ----------| ---------- | ----- |
|     1     |      77023 |  21.4 |
|     2     |      37338 |  10.4 |
|     4     |      22462 |   6.2 |
|     8     |      13430 |   3.7 |
|    10     |      10698 |   3.0 |
|    16     |       7077 |   2.0 |
|    24     |       4888 |   1.4 |
|    32     |       3713 |   1.0 |
|    64     |       5641 |   1.6 |

This was run on a machine used by others, so timings have likely been perturbed
by other loads on the machine.

Note that going from 32 to
64 processors significantly hurts performance, although in this case this
experiment revealed that the computer, despite reporting 128 cores appears to
either have just 32 processors, users might be limited to 32 cores... Or maybe
the code really cannot scale beyond 32 threads. Parallel computing can be
mystifying.

Gaphically this looks like
[spread_speedup](Images/spread_speedup.png)

This is typical of many OMP applications where the addition of more
threads will provide diminishing benefit as threads are added to the
task.  Not however that in this
case it does not exhibit significant worsening of the trend as threads
are added.  Different programmes have different parallel
characteristics.


# Summary
LillyMol has a variety of tools designed for local parallel, cluster
parallel and thread level parallel processing. While the local parallel
and cluster parallel uses are obvious the most efficient, the thread
parallel tools help with tasks that do not have obviously parallel
implementations.
