#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(readlink -e $0))))
fi

config=${LILLYMOL_HOME}/data/minor_changes.textproto

$LILLYMOL_HOME/bin/$(uname)/minor_changes -C ${config} "$@"
