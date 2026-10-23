#!/usr/bin/env bash

# Version 1 of a chirality fingerprint.
# Crude
# Pipes the input to cip_labeler which places isotopes on the marked chiral centres.
# Pipes that output to iwecfp with start atoms as the marked atoms.

argv=("$@")
input=${argv[-1]}
unset 'argv[-1]'

# echo "Input file is ${input}"
# echo "Other args"
# echo $argv

query_file="${LILLYMOL_HOME}/data/queries/chirality_fingerprint/iso89.qry"

if [[ ! -s "${query_file}" ]] ; then
  echo "Cannot find query file ${query_file}" >&1
  exit 1
fi

# If the -f option is passed, we need to add it to cip_labeler

dash_f=''
for opt in "${argv[@]}" ; do
  if [[ "${opt}" == "-f" ]] ; then
    dash_f='-f'
    break
  fi
done

cip_labeler.sh ${dash_f} ${input} | iwecfp.sh ${argv[@]} -q PROTO:${query_file} -
