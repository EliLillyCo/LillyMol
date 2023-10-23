#!/bin/sh

# Setup queries for jwcats

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi

charges="${LILLYMOL_HOME}/data/queries/charges/queries"
hbonds="${LILLYMOL_HOME}/data/queries/hbonds"

exe=${LILLYMOL_HOME}/bin/$(uname)/jwcats
exec ${exe} -N noremove -N F:${charges} -H noremove -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label "$@"

