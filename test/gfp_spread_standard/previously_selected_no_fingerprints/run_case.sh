#!/bin/bash

# gfp_spread_standard must fail, loudly, when the file given to -A yields no
# fingerprints. Before the fix it reported success and silently ignored -A,
# leaving every initial distance at 1.0 - indistinguishable from -A never
# having been specified - and produced a full spread anyway.
#
# The harness passes: $1 executable, $2 in dir, $3 out dir, then the args from
# tests.json ($4 here is the shared data directory).
#
# Everything echoed is deliberately free of file paths so that the golden
# output is portable between checkouts.

exe=$1
datadir=$4

if [[ -z ${exe} || -z ${datadir} ]] ; then
  echo 'usage: run_case.sh <exe> <indir> <outdir> <datadir>'
  exit 1
fi

for bad in empty.gfp not_tdt.gfp ; do
  "${exe}" -h 1 -A "${datadir}/${bad}" "${datadir}/pool.gfp" > spread.out 2> spread.err
  rc=$?

  if [[ ${rc} -eq 0 ]] ; then
    echo "${bad} exit_status zero"
  else
    echo "${bad} exit_status nonzero"
  fi

  echo "${bad} records_written $(grep -c '^|$' spread.out)"
  echo "${bad} diagnostic_lines $(grep -c 'no fingerprints read' spread.err)"
done

exit 0
