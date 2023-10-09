#!/bin/bash

if [[ ! -v LILLYMOL_HOME ]]
then
  echo 'Must set shell variable LILLYMOL_HOME' >&2
  exit 1
fi

charges="${LILLYMOL_HOME}/contrib/data/queries/charges/queries"
hbonds="${LILLYMOL_HOME}/contrib/data/queries/hbonds"
# File containing reduced_graph::GraphReduction testprotos describing
# a series of reductions.
reductions=${LILLYMOL_HOME}/contrib/data/reduced_graph/reductions

# If you want to use reductions based on LillyMol donor acceptor assignments, use
# reductions=${LILLYMOL_HOME}/contrib/data/reduced_graph/reductions.donor_acceptor
# Those queries use the isotopic labels applied by the donor acceptor.

exec ${LILLYMOL_HOME}/bin/Linux/reduced_graph -N F:${charges} \
          -R ${reductions} \
          -H noremove -H a=F:${hbonds}/acceptor -H d=${hbonds}/donor.qry -H label \
          -N noremove -N F:${charges} \
          "$@"
