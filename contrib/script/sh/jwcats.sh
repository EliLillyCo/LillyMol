#!/bin/sh

# Setup queries for jwcats

if [ ! -v LILLYMOL_HOME ] ; then
  echo "Must set shell variable LILLYMOL_HOME" >&2
  exit 1
fi

charges="${LILLYMOL_HOME}/contrib/data/queries/charges/queries"
hbonds="${LILLYMOL_HOME}/contrib/data/queries/hbonds"

exec ${LILLYMOL_HOME}/bin/Linux/jwcats -N noremove -N F:${charges} -H noremove -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label "$@"

