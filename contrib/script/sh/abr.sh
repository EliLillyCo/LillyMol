#!/bin/bash

# set up queries for abraham queries.

if [ ! -v LILLYMOL_HOME ] ; then
  echo "Must set shell variable LILLYMOL_HOME" >&2
  exit 1
fi

ABR_LIB=${LILLYMOL_HOME}/data/abraham

F=${ABR_LIB}/Abraham
P=${ABR_LIB}/Alpha2H
exec ${LILLYMOL_HOME}/bin/Linux/abraham -E autocreate -l -F ${F} -P ${P} -g all  "$@"

