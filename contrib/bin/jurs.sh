#!/usr/bin/env bash

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME="$(dirname $(realpath $0))/../../"
fi

executable="${LILLYMOL_HOME}/bin/$(uname)/jwsadb"

if [[ ! -x ${executable} ]] ; then
  executable='jwsadb'
fi

DATA_DIR=${LILLYMOL_HOME}/data
QUERIES=${DATA_DIR}/queries

exec $executable -N noremove -N F:${QUERIES}/charges/queries \
         -H noremove -H a=F:${QUERIES}/hbonds/acceptor \
         -H d=${QUERIES}/hbonds/donor.qry -H label \
         -F ${DATA_DIR}/wildman_crippen.dat \
         -A I -A C -h -j -e -s -a -V molvol -i ignore_bad_chiral "$@"

