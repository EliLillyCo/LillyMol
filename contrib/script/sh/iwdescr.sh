#!/bin/bash

if [ ! -v LILLYMOL_HOME ] ; then
  echo "Must set shell variable LILLYMOL_HOME" >&2
  exit 1
fi

charges="${LILLYMOL_HOME}/contrib/data/queries/charges/queries"
hbonds="${LILLYMOL_HOME}/contrib/data/queries/hbonds"
ranges=${LILLYMOL_HOME}/contrib/data/chembl.ranges

exec ${LILLYMOL_HOME}/bin/Linux/iwdescr -N F:${charges} \
        -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label \
        -B ranges=${ranges} \
        -l -g all -u 0 -b 5 -O all -B quiet -E autocreate "$@"

