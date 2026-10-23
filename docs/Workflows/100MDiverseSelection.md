# Diverse selection from 100M molecules.

## Objective
The problem involved selecting several thousand diverse molecules
from a collection of over 100M.

## Strategies
Generally I would proposed finding some means of prioritising the molecules, even
if it were just by number of heavy atoms. But in this case, there did not appear
to be any reasonable means of prioritising the molecules.

Unfortunately, none of the LillyMol tools can process 100M+ molecules in any
reasonable time - or memory.

Therefore approximate solutions are needed. 

Agree that any subset containing ~8 thousand molecules is very likely an inadequate
representation of 100M+ molecules, even if the 100M are combinatorially derived - 
containing several different libraries.

Generally gfp_spread is the preferred tool for selecting diverse molecules. But we
also observe that spread, correctly, identifies the most unusual molecules in any
collection. In this case, that was not an objective, so instead, a leader based
approach was adopted.

## One Solution
Given our inability to perform leader clustering on the whole set, subdivide the
set into 1M molecule chunks with iwsplit
```
iwsplit -suffix smi -n 100000 big.smi
```
generates a bunch of `iwsplit*.smi` files, each containing 1M molecules. For
each chunk, generate fingerprints and run leader. The -s option to leader
is needed to let it know how many fingerprints are in the input. Also
for each chunk, select only the first 8000 clusters.

Do this 16 way parallel on the current machine.
```
dopattern.sh -parallel 16 -suffix smi 'gfp_make.sh iwsplit%.smi | \
   gfp_leader_standard -n 8000 -h 1 -t 0.25 -s 1000000 - > %.ldr'
```
Note that we run each leader invocation with just one thread. Multi-threading
is the enemy of throughput, and this is all about throughput.

The choice of 0.25 for the leader radius was empirically determined.
Given 1M molecules per file, what is a threshold that allows sampling
of almost all of the molecules in that file. If the radius is too small,
then selecting 8K can likely be sampled by examining only a small
subset of the file. In this case, for these molecules, a radius of
0.25 generally allowed sampling across all of the file.

Selecting 8k molecules from 1M fingerprints with leader is relatively fast, 5 minutes or so
on a modern CPU (2025). Each of these tasks is fully independent, so can be
processed in parallel. In this particular instance, the tasks were done
on a 16 core server, but they could also have been done via a cluster.

Total run time on the 16 core server was under two hours.

This yields 100+ leader output files, each containing the first 8k
clusters from that file. The molecules within these files are internally
diverse, with no pair of molecules being closer than 0.25.

Convert these to smiles files containing only the cluster centre.
In this case it is safer to use the number of splits for dopattern,
because smiles files are being produced.

```
dopattern.sh -o 150 'nplotnn -L def -n 0 iwsplit%.ldr > iwsplit%.ldr.smi'
```

While each file is internally diverse, we can use leader to select
a final set of 8k molecules from it.

```
gfp_make.sh iwsplit*.ldr.smi > final.gfp
gfp_leader_standard -h 8 -n 8000 -t 0.25 final.gfp > final.ldr
```
This time, we use multi-threading to perform the final 8k
candidates, throughput is no longer a concern.

## Analysis and Optimisations
While the final selections will be internally diverse, it is
unclear by how much the final outcome would be improved if
it were possible to run leader on the starting 100M molecules.

