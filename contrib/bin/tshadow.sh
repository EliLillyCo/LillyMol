#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME="$(dirname $(realpath $0))/../../"
fi

executable="${LILLYMOL_HOME}/bin/Linux/tshadow"

if [[ ! -x ${executable} ]] ; then
  executable='tshadow'
fi

exec ${executable} -L -G "$@"
