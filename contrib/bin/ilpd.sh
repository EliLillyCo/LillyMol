#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  echo "Must set $LILLYMOL_HOME" >&2
  exit 1
fi

datadir="${LILLYMOL_HOME}/contrib/data/ILPD"
here=$(dirname $0)

config="${datadir}/config.textproto"

if [[ ! -s "${config}" ]] ; then
  echo "Configuration file ${config} not missing or empty"
  exit 1
fi

extra_charge="${datadir}/extra_charge.qry"

exec ${LILLYMOL_HOME}/bin/$(uname)/ilpd_reduced_graph -C "${config}" -N F:${LILLYMOL_HOME}/data/queries/charges/queries \
        -N minsep=1 -N ${extra_charge} "$@"
