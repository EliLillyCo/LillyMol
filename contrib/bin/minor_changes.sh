#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(realpath $0))))
fi

config=${LILLYMOL_HOME}/data/minor_changes/minor_changes.textproto
fragments=${LILLYMOL_HOME}/data/minor_changes/fragments.textproto

$LILLYMOL_HOME/bin/$(uname)/minor_changes -P UST:AY -F ${fragments} -C ${config} "$@"
