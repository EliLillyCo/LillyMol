#!/bin/bash

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi

charges="${LILLYMOL_HOME}/data/queries/charges/queries"
hbonds="${LILLYMOL_HOME}/data/queries/hbonds"
ranges=${LILLYMOL_HOME}/data/chembl.ranges

exec ${LILLYMOL_HOME}/bin/Linux/iwdescr -N F:${charges} \
        -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label \
        -B ranges=${ranges} \
        -l -g all -u 0 -b 5 -O all -B quiet -E autocreate "$@"

