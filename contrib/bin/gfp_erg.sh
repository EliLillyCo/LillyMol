#!/bin/bash

# Set up queries for gfp_erg
here=$(dirname $(readlink -e $0))

if [[ -v LILLYMOL_HOME ]] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname ${here}))
fi

ERG_QUERIES=${LILLYMOL_HOME}/data/ErG

if [ ! -d ${ERG_QUERIES} ] ; then
  echo "Where is ${ERG_QUERIES}" >&2
  exit 1
fi

CHARGES="-N F:${ERG_QUERIES}/charges/queries"
ACCEPTOR="-H a=F:${ERG_QUERIES}/hbonds/acceptor"
DONOR="-H d=${ERG_QUERIES}/hbonds/donor_nik.qry"
LABEL="-H label"
ENDCAPS="-q F:${ERG_QUERIES}/endcap/endcaps"
FUZZY="-F 0.3"
SET="-s"
DIST="-d 10"
ABSTRACT="-m"
MAX_FF="-X 20"

exe=${LILLYMOL_HOME}/bin/$(uname)/gfp_erg
exec ${exe} -n ${AROMATICITY} ${CHARGES} ${ACCEPTOR} ${DONOR} ${LABEL} ${ENDCAPS} ${FUZZY} "$@"
