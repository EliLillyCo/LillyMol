#!/usr/bin/env bash

# set up queries for abraham queries.

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi

executable="${LILLYMOL_HOME}/bin/$(uname)/dbf"
if [ ! -s "${executable}" ]; then
    echo "Cannot access executable '$executable'" >&2
    exit 1
fi

QUERIES="${LILLYMOL_HOME}/data/queries"

exec $executable \
        -N F:${QUERIES}/charges/queries \
        -H a=F:${QUERIES}/hbonds/acceptor \
        -H d=${QUERIES}/hbonds/donor.qry \
        -H label -s '[*+]' -s '[*-]' -s '[1*,2*]' -s '[2*,3*]' \
        -u 0 -p 2d -p 3d -r -A D -A I -A C \
        -i ignore_bad_chiral -i sdf "$@"
