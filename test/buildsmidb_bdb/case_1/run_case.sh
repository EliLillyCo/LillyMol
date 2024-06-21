#!/usr/bin/env bash

me=$0
me=$(readlink -e $0)

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME=$(dirname $(dirname $(dirname $(dirname ${me}))))
fi

if [[ ! -v BUILD_DIR ]] ; then
  BUILD_DIR=$(uname)
fi

build=${LILLYMOL_HOME}/bin/${BUILD_DIR}/buildsmidb_bdb
lookup=${LILLYMOL_HOME}/bin/${BUILD_DIR}/in_database_bdb

if [[ ! -x ${build} || ! -x ${lookup} ]] ; then
  echo "${build} and/or ${lookup} not available" >&2
  exit 1
fi

insmi="$(dirname ${me})/in/in.smi"

if [[ ! -s ${insmi} ]] ; then
  echo "MIssing input file #{insmi}" >&2
  exit 1
fi

dbname="/tmp/buildsmidb$$.bdb"

cmd="${build} -d ${dbname} -c -l ${insmi}"
${cmd}

if [[ $? -ne 0 ]] ; then
  echo "${cmd} failed" >&2
  exit 1
fi

if [[ ! -s ${dbname} ]] ; then
  echo "Did not build database ${dbname}" >&2
  exit 1
fi

found="/tmp/found$$"
not_in_db="/tmp/notfound$$"
cmd="${lookup} -d ${dbname} -c -l -F ${found} -U ${not_in_db} ${insmi}"
${cmd}

if [[ $? -ne 0 ]] ; then
  echo "${cmd} failed" >&2
  exit 1
fi

if [[ -s "${not_in_db}.smi" ]] ; then
  echo "Some molecules not found in db" >&2
  exit 1
fi

if [[ ! -s "${found}.smi" ]] ; then
  echo "DId not find any molecules in db" >&1
  exit 1
fi

unlink ${dbname}
unlink "${found}.smi"
unlink "${not_in_db}.smi"
exit 0
