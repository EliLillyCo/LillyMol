#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  here=$(basename $0)
  export LILLYMOL_HOME=${here}/../..
fi

charges=${LILLYMOL_HOME}/data/queries/charges/queries
hbonds=${LILLYMOL_HOME}/data/queries/hbonds
donor_acceptor="-H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label"
hbondpatterns="${LILLYMOL_HOME}/data/queries/hbondpatterns"
exec "${LILLYMOL_HOME}/bin/$(uname)/tnass" -a -a -N F:${charges} ${donor_acceptor} -q F:${hbondpatterns}/nass $@
