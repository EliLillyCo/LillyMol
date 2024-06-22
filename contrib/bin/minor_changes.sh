#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  me=$0
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname ${me})))
fi

config=${LILLYMOL_HOME}/data/minor_changes/minor_changes.textproto

$LILLYMOL_HOME/bin/$(uname)/minor_changes -C ${config} "$@"
