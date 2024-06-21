# GFP Tools

For many common tasks associated with molecular similarity there will be a
`gfp_*` tool. Sometimes there may be a serial version of that tool as well
as one or more parallel versions. This document should help you decide which
tool is most appropriate for a given circumstance.

Unless you are dealing with > tens of thousands of fingerprints, the serial
versions are likely to be adequate. For some comparion functions the
serial version may still be the best choice if dealing with > 1M molecules.

## Standard
Over time, we determined that a combination of 4 fingerprints generated
very good concordance with human perceptions of similarity. This is generated
via
```
gfp_make.sh -STD file.smi > file.gfp
```
and is also the default. Whereas the GFP framework is designed to work with
arbitrary combinations of fingerprints, we built some tools that only consume
the 'standard' set, in order to gain some computational advantages. Many of
the 'standard' tools are also parallelised.

## Intra-Collection similarity
To find the pair-wise similarities within a collection of molecules use
```
gfp_nearneighbours_single_file
gfp_nearneighbours_single_file_tbb
```
The first is a serial version, while the second is parallelised
using Intel's Threading Building Blocks (TBB)[https://github.com/oneapi-src/oneTBB].

See [gfp_nearneighbours_single_file](gfp_nearneighbours_single_file.md)

## Inter-Library similarity
There are several tools for doing this. Generally we will be comparing a
small set of molecules against a possibly much larger set of molecules.
In these scenarios we refer to the two sets as the 'needles' and the 'haystack'.

Generally the preferred tools are
```
gfp_lnearneighbours
gfp_lnearneighbours_standard
```
These work by reading the needles into memory - via the -p option. The
haystack is then read one at a time and compared with the needles. The neighbours
are accumulated. Once the haystack has been read, the neighbours are written.
Often the haystack will be read from stdin.

The 'standard' version runs parallel, so is best used with larger sets of
needles. But note that you might ultimately get better throughput by using
smaller sets of needles, spread across more parallel jobs.

## Diverse Selection (spread)
This is a greedy max/min distance selector, that will result in a set of
molecules ranked from most to least diverse. Frequently used for selecting a
maximally diverse subset of molecules from a larger collection.

Tools are
```
gfp_spread
gfp_spread_standard
gfp_spread_omp
```
If using the standard fingerprints, gfp_spread_standard is by far the most
performant version, processing 50k fingerprints in 8 seconds using 8 cores.
gfp_spread_omp performs the same task in 32 seconds.

The default version, gfp_spread, performs the same task in 256 seconds,
indicating that the parellel versions are achieving good speedups.

The variant gfp_spread_buckets is a variant that enforces groups of
molecules to be selected. TODO:ianwatson document this.

## Leader - Sphere Exclusion
This is commonly used for selecting a diverse yet desirable set of
molecules. Molecules are ordered by desirability, and leader run on
that set of fingerprints.

Tools are
```
gfp_leader
gfp_leader_standard
gfp_leader_tbb
```

Again, gfp_leader runs serial, and can consume any kind of fingerprint, or
fingerprint combination. Given 50 random CHembl molecules clustered at a
0.25 radius, that task takes 126 seconds. The parallel version, gfp_leader_tbb
performs the same task in 22 seconds, while if you can use the standard fingerprints,
the task completes in 6 seconds.

