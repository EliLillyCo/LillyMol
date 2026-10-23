#!/usr/bin/env bash

here=$(dirname $0)

charge_assigner=${LILLYMOL_HOME}/data/queries/charges/queries
exec ${LILLYMOL_HOME}/bin/$(uname)/iwfp -a -a -e -p chgfp -O 7 -N F:${charge_assigner} "$@"
