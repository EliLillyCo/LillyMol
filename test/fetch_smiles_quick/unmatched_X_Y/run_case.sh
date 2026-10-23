#!/bin/bash
# Args from run_all_test.rb: executable, indir, outdir, then TestCase args.
# The -Y file is written by iterating a hash, so its order is not guaranteed.
# Sort it so the comparison is stable across platforms.
exe="$1"; datadir="$4"
"${exe}" -X notInIds -Y notInHaystack "${datadir}/ids_missing.txt" "${datadir}/haystack.smi" || exit $?
sort -o notInHaystack notInHaystack
