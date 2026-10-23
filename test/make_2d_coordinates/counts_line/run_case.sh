#!/bin/bash
# Args from run_all_test.rb: executable, indir, outdir, then TestCase args.
# LillyMol writes an abbreviated molfile counts line by default, which several
# other toolkits refuse to read. Coordinates that only LillyMol can read are no
# use to a depiction, so make_2d_coordinates asks for full V2000 records, with
# -t to get the terse form back. Only the counts line of each record is compared
# - it is platform independent, unlike the coordinates.
exe="$1"; datadir="$4"

"${exe}" -S full "${datadir}/input.smi" || exit $?
"${exe}" -t -S terse "${datadir}/input.smi" || exit $?

# The counts line is the 4th line of each record.
{
  echo 'default, full V2000:'
  awk '/^\$\$\$\$/{n=0;next} {n++} n==4' full.sdf
  echo 'with -t, terse:'
  awk '/^\$\$\$\$/{n=0;next} {n++} n==4' terse.sdf
} > counts_lines
