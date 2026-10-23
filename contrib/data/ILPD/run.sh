#!/usr/bin/env bash

here=$(dirname $0)

config="${here}/config.textproto"

if [[ ! -s "${config}" ]] ; then
  echo "Config file #{config} missing" >&2
  exit 1
fi

extra_charge="${here}/extra_charge.qry"

if [[ ! -s "${extra_charge}" ]] ; then
  echo "Extra charge query missing" >&2
  exit 1
fi

charge_assigner="${LILLYMOL_HOME}/data/queries/charges"

exec ipld_reduced_graph -C "${config}" -N ${charge_assigner}/queries \
        -N miinsep=1 -N ${extra_charge} "$@"

