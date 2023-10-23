#!/bin/bash

# set up queries for abraham queries.

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi

ABR_LIB=${LILLYMOL_HOME}/data/abraham
if [[ ! -d ${ABR_LIB} ]] ; then
  echo "$0 where are my queries ${ABR_LIB}" >&2
  exit 1
fi


F=${ABR_LIB}/Abraham
P=${ABR_LIB}/Alpha2H
exec ${LILLYMOL_HOME}/bin/Linux/abraham -E autocreate -l -F ${F} -P ${P} -g all  "$@"

