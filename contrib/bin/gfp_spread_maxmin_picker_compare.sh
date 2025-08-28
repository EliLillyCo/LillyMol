#!/usr/bin/env bash

# Back in 2017 Greg Landrum published an update introducing work Roger Sayle had done
# to improve the performance of the MaxMinPicker in RDKit.
# That functionality has been in LillyMol since the mid 1990's - it was one of the first
# projects I worked on at Lilly.

# The timings claimed by the new RDKit functionality sounds impressive and I was curious
# to see how LillyMol might compare.

# Use the last test case they cite.
# 247477 as the compounds that we already possess.
# 12386 compounds from which we wish to select

# Note that this benchmark uses a LillyMol default fingerprint that consists of four
# components.
# 1. a 2048 bit linear, path based fingerprint.
# 2,3. two 196 bit dictionary based fingerprints.
# 4. 8 molecular properties that contribute via ratios.

# So again, this makes comparisons even less valid since LillyMol is doing 
# more work - evaluating the properties tends to be somewhat slow in comparison
# to fingerprints.

# We don't have the exact molecules so in each case we randomly select the appropriate
# number from Chembl - which must be the first argument to this script

if [[ "$#" == 0 ]] ; then
  echo "Must specify path name of Chembl - or some other large collection of molecules"
  exit 1
fi

haystack=$1

tmpdir='/tmp'

random_records -n 247477 ${haystack} > "${tmpdir}/247477.smi"
random_records -n 12386 ${haystack} > "${tmpdir}/12386.smi"

echo "Generating fingerprints 247k"
time gfp_make.sh ${tmpdir}/247477.smi > ${tmpdir}/247477.gfp
echo "Generating fingerprints 12k"
time gfp_make.sh ${tmpdir}/12386.smi > ${tmpdir}/12386.gfp

echo "Select 1 molecule no parallelism, cmp 8.5"
time gfp_spread_standard -h 1 -n 1 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/1.spr

echo "Select 1 molecule 2 threads"
time gfp_spread_standard -h 2 -n 1 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/1h2.spr

echo "Select 1 molecule 8 threads"
time gfp_spread_standard -h 8 -n 1 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/1h8.spr

echo "Select 2 molecules no parallelism, cmp 133"
time gfp_spread_standard -h 1 -n 2 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/2h1.spr

echo "Select 2 molecules 2 threads"
time gfp_spread_standard -h 2 -n 2 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/2h2.spr

echo "Select 2 molecules 8 threads"
time gfp_spread_standard -h 8 -n 2 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/2h8.spr

echo "Select 3 molecules no parallelism, cmp 218"
time gfp_spread_standard -h 1 -n 3 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/3h1.spr

echo "Select 3 molecules 2 threads"
time gfp_spread_standard -h 2 -n 3 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/3h2.spr

echo "Select 3 molecules 8 threads"
time gfp_spread_standard -h 8 -n 3 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/3h8.spr

echo "Select 1000 molecules no parallelism, cmp 1187"
time gfp_spread_standard -h 1 -n 1000 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/1000h1.spr

echo "Select 1000 molecules 2 threads"
time gfp_spread_standard -h 2 -n 1000 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/1000h2.spr

echo "Select 1000 molecules 8 threads"
time gfp_spread_standard -h 8 -n 1000 -A ${tmpdir}/247477.gfp ${tmpdir}/12386.gfp > /tmp/1000h8.spr

# On a consumer desktop from 2021
# 12th Gen Intel(R) Core(TM) i7-12700K
# the following timings are observed.

#    Generating fingerprints 247k
#    read mol smi eof
#    
#    real    0m17.975s
#    user    0m36.904s
#    sys     0m0.271s
#    Generating fingerprints 12k
#    read mol smi eof
#    
#    real    0m0.962s
#    user    0m1.882s
#    sys     0m0.026s
#    Select 1 molecule no parallelism, cmp 8.5
#    
#    real    0m43.782s
#    user    0m43.708s
#    sys     0m0.055s
#    Select 1 molecule 2 threads
#    
#    real    0m23.185s
#    user    0m44.724s
#    sys     0m0.640s
#    Select 1 molecule 8 threads
#    
#    real    0m8.484s
#    user    0m52.409s
#    sys     0m4.190s
#    Select 2 molecules no parallelism, cmp 133
#    
#    real    0m43.720s
#    user    0m43.639s
#    sys     0m0.071s
#    Select 2 molecules 2 threads
#    
#    real    0m23.244s
#    user    0m44.791s
#    sys     0m0.710s
#    Select 2 molecules 8 threads
#    
#    real    0m8.604s
#    user    0m52.997s
#    sys     0m4.248s
#    Select 3 molecules no parallelism, cmp 218
#    
#    real    0m43.477s
#    user    0m43.394s
#    sys     0m0.064s
#    Select 3 molecules 2 threads
#    
#    real    0m23.149s
#    user    0m44.678s
#    sys     0m0.651s
#    Select 3 molecules 8 threads
#    
#    real    0m8.505s
#    user    0m52.497s
#    sys     0m4.262s
#    Select 1000 molecules no parallelism, cmp 1187
#    
#    real    0m43.852s
#    user    0m43.767s
#    sys     0m0.067s
#    Select 1000 molecules 2 threads
#    
#    real    0m23.319s
#    user    0m44.960s
#    sys     0m0.669s
#    Select 1000 molecules 8 threads
#    
#    real    0m8.672s
#    user    0m53.483s
#    sys     0m4.078s

# While comparisons are impossible due to different hardware, in most cases we
# see significantly faster times with gfp_spread_standard, and that the multi
# threading option seems to be effective - on the last case we see time dropping
# from 43.9 seconds to 8.7 by using 8 cores.

# On CPU's with more threads, using 16 cores can show further improvement.

# Roger's implementation of selecting the first molecule shows a very effective
# optimisation not present in gfp_spread_standard.

# Note that selecting 1 or 1000 molecules using gfp_spread_standard takes about
# the same time. The times are different, but the difference between selecting
# one molecule and 1000 is about 0.2 seconds on this machine.

# The last computation, selecting 1000 molecules with 8 cores is done in 8.7
# seconds. This compares with 1187 from RDKit. I am not sure why there is such
# large difference, a difference that I don't think can be explained by hardware
# differences. Even without parallelism, the 43.9 seconds observed here is
# much better than the 1187 seen with RDKit.

# On a machine with more cores
# Intel(R) Xeon(R) CPU Max 9462
#
# selecting 1000 molecules with  8 threads takes 9.7 seconds.
# selecting 1000 molecules with 16 threads takes 6.5 seconds.
# selecting 1000 molecules with 24 threads takes 5.8 seconds.

# Clearly into the realm of significantly diminishing returns.
