#!/bin/bash

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(readlink -e $0))))
fi

exec ${LILLYMOL_HOME}/bin/Linux/hydrophobic_sections -E autocreate -G def -L def -i smi "$@"
