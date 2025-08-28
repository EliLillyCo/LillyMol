#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME="${0}/../.."
fi

config="${LILLYMOL_HOME}/data/alogp/aclogp.textproto"

exec ${LILLYMOL_HOME}/bin/Linux/alogp -C "${config}" "$@"
