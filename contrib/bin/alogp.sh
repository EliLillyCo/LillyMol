#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME || ! -d "${LILLYMOL_HOME}"  ]]; then
  echo "LILLYMOL_HOME not set" >&2
  exit 1
fi

QUERIES="${LILLYMOL_HOME}/data/queries"

charges="-N F:${QUERIES}/charges/queries"

exec $LILLYMOL_HOME/bin/$(uname)/alogp ${charges} -Y alcacid -Y RDKIT.N+ "$@"
