#!/bin/bash

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(readlink -e $0))))
fi

charges="${LILLYMOL_HOME}/data/queries/charges/queries"
hbonds="${LILLYMOL_HOME}/data/queries/hbonds"
# File containing reduced_graph::GraphReduction textprotos describing
# a series of reductions.
reductions=${LILLYMOL_HOME}/data/reduced_graph/reductions

# If you want to use reductions based on LillyMol donor acceptor assignments, use
# reductions=${LILLYMOL_HOME}/data/reduced_graph/reductions.donor_acceptor
# Those queries use the isotopic labels applied by the donor acceptor.

exec ${LILLYMOL_HOME}/bin/Linux/reduced_graph -N F:${charges} \
          -R ${reductions} \
          -H noremove -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label \
          -N noremove -N F:${charges} \
          "$@"
