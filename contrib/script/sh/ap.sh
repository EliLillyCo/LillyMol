#!/bin/bash

# set up queries for Abraham Platts queries.

if [ ! -v LILLYMOL_HOME ] ; then
  echo "Must set shell variable LILLYMOL_HOME" >&2
  exit 1
fi

CONSTANTS=${LILLYMOL_HOME}/data/abraham_hbond_constants

exec ${LILLYMOL_HOME}/bin/Linux/substituent_model -E autocreate -l \
        -q s:${CONSTANTS} -u -M ap_sa= -M ap_sb= "$@"
