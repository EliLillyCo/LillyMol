#!/bin/bash

# Setup queries for jwcats

if [[ -v LILLYMOL_HOME ]] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $0)))
fi

charges="${LILLYMOL_HOME}/data/queries/charges/"
hbonds="${LILLYMOL_HOME}/data/queries/hbonds"

exec ${LILLYMOL_HOME}/bin/Linux/jwdist \
     -N athpos -N F:${charges}/queries \
     -q QUERY:F:${charges}/positive \
     -q QUERY:F:${charges}/negative \
     -q SMARTS:'[!#6&!#1]' -q SMARTS:'[c,n,o]' \
     -q SMARTS:'[Cl,Br,I,F]' -q SMARTS:'[C,N,O,S]=[C,N,O,S]' \
     -q SMARTS:'[#1][N,O,P,S]' -q SMARTS:'[!#1]' \
     -r 10 -S - -A I -A C -i ignore_bad_chiral -i sdf "$@"
