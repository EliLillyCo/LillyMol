#!/usr/bin/env bash

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(realpath $0))))
fi

exec ${LILLYMOL_HOME}/bin/$(uname)/hydrophobic_sections -E autocreate -G def -L def "$@"
