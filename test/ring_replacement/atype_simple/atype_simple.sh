#!/usr/bin/env bash

exe="$1"
indir="$2"
outdir="$3"

bindir=$(dirname "$exe")
ring_extraction="${bindir}/ring_extraction"
ring_replacement="${bindir}/ring_replacement"

atype='-P UST:Z'

${ring_extraction} ${atype} -S RINGS ${indir}/library.smi

rings6a='RINGS_6a.smi'

if [[ ! -s ${rings6a} ]] ; then
  echo "Did not create RINGS_6a.smi" >&2
  exit 1
fi

${ring_replacement} -R ${rings6a} ${atype} ${indir}/in.smi
