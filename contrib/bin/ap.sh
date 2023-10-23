#!/bin/bash

# set up queries for Abraham Platts queries.

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi
echo "LILLYMOL_HOME in ${LILLYMOL_HOME}"

CONSTANTS=${LILLYMOL_HOME}/data/abraham_hbond_constants
ls ${CONSTANTS}

exec ${LILLYMOL_HOME}/bin/Linux/substituent_model -E autocreate -l \
        -q S:${CONSTANTS} -u -M ap_sa= -M ap_sb= "$@"
