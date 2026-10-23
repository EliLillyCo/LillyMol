#!/usr/bin/env bash

here=$(dirname $(realpath $0))

if [[ -v LILLYMOL_HOME ]] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(realpath $0)))
fi

exec ruby ${here}/dopattern.rb "$@"
