#!/bin/bash

here=$(dirname $(readlink -e $0))

if [[ -v LILLYMOL_HOME ]] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(readlink -e $0)))
fi

exec ruby ${here}/dopattern.rb "$@"
