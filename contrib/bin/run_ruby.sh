#!/usr/bin/env bash
here=$(dirname $(realpath $0))

if [[ ! -v LILLYMOL_HOME ]] ; then
  export LILLYMOL_HOME=$(dirname $(dirname $here))
  export PATH=${here}:${PATH}
fi

me=$(basename $0)
exec ruby ${here}/${me/.sh/.rb} "$@"
