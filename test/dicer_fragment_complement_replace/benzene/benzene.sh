#!/usr/bin/env bash

exe="$1"
indir="$2"
outdir="$3"

bindir=$(dirname "$exe")
dicer="${bindir}/dicer"

dicer_options='-B brcb -C iso=1'

serialized='start.diced.bin'
${dicer} ${dicer_options} -B serialized_proto -S ${serialized} ${indir}/start.smi
if [[ ! -s ${serialized} ]] ; then
  echo "Did not create ${serialized}" >&2
  exit 1
fi


dicer_to_complemenent_db="${bindir}/dicer_to_complemenent_db"

bdb="test_benzene.bdb"
${dicer_to_complemenent_db} -d ${bdb} ${serialized}


if [[ ! -s ${bdb} ]] ; then
  echo "BDB not created" >&2
  exit 1
fi


# dice the input to generate textproto
diced_textproto='benzene.diced.textproto'
${dicer} ${dicer_options} -B proto ${indir}/start.smi > ${diced_textproto}

if [[ ! -s ${diced_textproto} ]] ; then
  echo "Could not create diced textproto" >&2
  exit 1
fi

dicer_fragment_complement_replace="${bindir}/dicer_fragment_complement_replace"
${dicer_fragment_complement_replace} -p PARENT -v -d ${bdb} ${diced_textproto}
