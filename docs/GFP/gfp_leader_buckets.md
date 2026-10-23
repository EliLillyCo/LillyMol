# gfp_leader_buckets

`gfp_leader_buckets` performs leader, or sphere exclusion, clustering on
fingerprints while trying to preserve representation from a set of user-defined
buckets.

This is useful when a collection has already been divided into classes, assays,
clusters, model-score bins, property bins, or other groups, and the selected
cluster centres should sample those groups as uniformly as possible.

The tool reads one GFP/TDT fingerprint file and writes leader-cluster TDT records
to standard output.

## Quick Start

Generate fingerprints and use the second token of the identifier as the bucket
label.

```shell
gfp_make.sh input.smi > input.gfp
gfp_leader_buckets -t 0.25 -B col=2 -n 500 input.gfp > selected.ldr
```

Use an external id-to-bucket file.

```shell
gfp_leader_buckets -t 0.25 -B buckets.txt -n 500 input.gfp > selected.ldr
```

Write a per-bucket summary.

```shell
gfp_leader_buckets -t 0.25 -j 8 -Y bucket_summary.csv -B buckets.txt input.gfp > selected.ldr
```

The selected leaders can be extracted from the leader output with `nplotnn`.

```shell
nplotnn -n 0 selected.ldr > leaders.smi
```

## Algorithm

The normal `gfp_leader` algorithm starts at the front of the input file, selects
the first available fingerprint as a leader, and removes all available
fingerprints within the threshold from further consideration.

`gfp_leader_buckets` keeps that same sphere exclusion step, but changes how the
next leader is chosen. Each candidate belongs to one bucket. At each selection
step, the tool examines bucket state and chooses the next leader from the bucket
that currently has the fewest leaders. Ties are resolved by choosing the bucket
with the fewest remaining available candidates. Within the chosen bucket, the
first still-available fingerprint in input order is selected.

This means the input order still matters within each bucket. Place the most
desirable molecules earlier in the input if the selection should prefer them.

Once a leader has been selected, all available fingerprints within `-t` of that
leader become cluster members and are no longer available as future leaders.
Bucket counts are updated after each cluster is formed.

## Options

```text
Usage: gfp_leader_buckets -t <dist> -B col=<n>|-B <file> [options] input.gfp

 -t <dist>      sphere exclusion distance threshold.
 -B col=<n>     bucket label is token <n> in the identifier/name field.
 -B <file>      two column id bucket mapping file; file=<fname> also accepted.
 -n <n>         maximum number of leaders to write.
 -j <n>         number of threads for distance scans, default 1.
 -r <n>         report progress every <n> leaders formed.
 -Y <fname>     write bucket summary; .csv suffix gives comma-separated output.
 -A <file>      previously selected fingerprints; candidates within threshold are discarded.
 -a <dist>      distance threshold for -A files; default is -t value.
 -I <tag>       identifier tag, default PCN<.
 -D ...         miscellaneous options, enter '-D help' for details.
 -F -P -W       standard GFP options, enter '-F help' for details.
 -v             verbose output.
```

### `-t <dist>`

Required. Specifies the sphere exclusion threshold. The value must be in the
range `[0,1]`.

No two selected leaders will be within this distance of each other, except for
any molecules that were supplied in a previously selected set via `-A`.

### `-B col=<n>`

Use token `<n>` in the identifier field as the bucket label. Token numbering is
1-based.

For example, with a fingerprint record containing

```text
PCN<CHEMBL12345 kinase>
```

`-B col=2` assigns this fingerprint to bucket `kinase`.

### `-B <file>`

Read bucket labels from a two-column text file.

```text
CHEMBL12345 kinase
CHEMBL67890 gpcr
CHEMBL24680 ion_channel
```

The first column is the fingerprint identifier and the second column is the
bucket label. Blank lines and lines starting with `#` are ignored. The form
`-B file=buckets.txt` is also accepted.

Every input fingerprint must have a bucket assignment when an external bucket
file is used. Missing assignments are fatal.

### `-n <n>`

Stop after writing at most `<n>` leaders. Without `-n`, clustering continues
until no available candidates remain.

### `-A <file>` and `-a <dist>`

`-A` supplies one or more files of previously selected fingerprints. Before
normal clustering begins, candidates within the specified threshold of any
previously selected fingerprint are discarded from the candidate pool.

By default, the previously selected threshold is the same as `-t`. Use `-a` to
specify a different threshold.

```shell
gfp_leader_buckets -t 0.25 -a 0.15 -A existing.gfp -B buckets.txt input.gfp > selected.ldr
```

Use `-a 0.0` when only exact fingerprint duplicates should be removed by the
previously selected set.

### `-j <n>`

Use `<n>` worker threads for distance scans. The clustering state is still
updated by the main thread, so the output is deterministic for a fixed input
file and command line.

Threading parallelizes the expensive search for neighbours within the threshold.
It does not change the leader selection rule.

For small remaining candidate lists, thread creation overhead can cost more than
the parallel scan saves. By default, `gfp_leader_buckets` runs the distance scan
serially unless at least 10000 candidates remain to be examined. Use
`-D minscan=<n>` to tune that threshold.


### `-r <n>`

Report progress to stderr every `<n>` leaders formed. The output TDT stream is
unchanged. This is useful for long timing runs where monitoring stdout size is
less convenient than watching a progress message.

### `-Y <fname>`

Write a bucket summary file. If `<fname>` ends in `.csv`, the summary is comma
separated. Otherwise it is space separated.

The columns are

```text
bucket total leaders cluster_members previously_selected available
```

where

```text
total = leaders + cluster_members + previously_selected + available
```

`previously_selected` counts candidates discarded because of an `-A` file.
`available` will usually be zero when clustering runs to completion, but may be
non-zero when `-n` stops the run early.

## Output

The output is a TDT stream in the usual leader-clustering form. Each cluster
contains the leader followed by the cluster members found within the threshold.
The bucket label is written with the `BCKT<...>` dataitem.

```text
$SMI<CCO>
PCN<mol1>
CLUSTER<0>
CSIZE<2>
BCKT<bucket_a>
$SMI<CCN>
PCN<mol2>
DIST<0.18>
BCKT<bucket_a>
|
```

Cluster members are written in input order. `DIST<...>` is written for cluster
members and is the distance from the leader.

## Notes

`gfp_leader_buckets` is intended for bucket-balanced selection. If no bucket
balancing is needed, use `gfp_leader`.

The algorithm can only balance buckets to the extent allowed by the chemistry
and the threshold. A small bucket may be exhausted quickly if most of its
molecules are close to selected leaders, or if many of them are discarded by a
previously selected set.

Very small thresholds tend to select more leaders and give the bucket balancing
rule more opportunities to act. Larger thresholds form larger clusters and may
exhaust some buckets early.
