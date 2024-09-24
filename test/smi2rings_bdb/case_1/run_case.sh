#!/usr/bin/env bash

me=$0
me=$(readlink -e $0)

if [[ ! -v LILLYMOL_HOME ]] ; then
  LILLYMOL_HOME=$(dirname $(dirname $(dirname $(dirname ${me}))))
fi

if [[ ! -v BUILD_DIR ]] ; then
  BUILD_DIR=$(uname)
fi

diff_tool='../../fileDiff.sh'

smi2rings=${LILLYMOL_HOME}/bin/${BUILD_DIR}/smi2rings_bdb

if [[ ! -x ${smi2rings} ]] ; then
  echo "${smi2rings} not available" >&2
  exit 1
fi

echo "Testing:  $smi2rings"

insmi="$(dirname ${me})/in/in.smi"

if [[ ! -s ${insmi} ]] ; then
  echo "Missing input file #{insmi}" >&2
  exit 1
fi

dbname="/tmp/smi2rings$$.bdb"

stderr='stderr'

common_options="-d ${dbname} -n -g all -j double -j spiro -z 8 -N add -c -v"
cmd="${smi2rings} -d STORE -j ring -j iso -j env=UST:achry ${common_options} ${insmi}"
${cmd} 2> ${stderr}

if [[ $? -ne 0 ]] ; then
  echo "${cmd} failed" >&2
  exit 1
fi

if [[ ! -s ${dbname} ]] ; then
  echo "Did not build database ${dbname}" >&2
  exit 1
fi

found="/tmp/found$$"

cmd="${smi2rings} -d LOOKUP -a ${common_options} -j ring ${insmi}"
${cmd} > ${found} 2> ${stderr}

if [[ $? -ne 0 ]] ; then
  echo "${cmd} failed" >&2
  exit 1
fi

if [[ ! -s "${found}" ]] ; then
  echo "Did not find any molecules in db" >&1
  exit 1
fi

${diff_tool} ${found} out/found_ring.smi
if [[ $? -eq 1 ]] ; then
  echo "case ring : TEST PASS"
else
  echo "case ring : TEST FAIL"
fi

cmd="${smi2rings} -d LOOKUP -a ${common_options} -j iso ${insmi}"
${cmd} > ${found} 2> ${stderr}

if [[ $? -ne 0 ]] ; then
  echo "${cmd} failed" >&2
  exit 1
fi

if [[ ! -s "${found}" ]] ; then
  echo "Did not find any molecules in db" >&1
  exit 1
fi

${diff_tool} ${found} out/found_iso.smi
if [[ $? -eq 1 ]] ; then
  echo "case iso : TEST PASS"
else
  echo "case iso : TEST FAIL"
fi

cmd="${smi2rings} -d LOOKUP -a ${common_options} -j env=UST:achry ${insmi}"
${cmd} > ${found} 2> ${stderr}

if [[ $? -ne 0 ]] ; then
  echo "${cmd} failed" >&2
  exit 1
fi

if [[ ! -s "${found}" ]] ; then
  echo "Did not find any molecules in db" >&1
  exit 1
fi

${diff_tool} ${found} out/found_env.smi
if [[ $? -eq 1 ]] ; then
  echo "case env : TEST PASS"
else
  echo "case env : TEST FAIL"
fi

unlink ${dbname}
unlink "${found}"
unlink "${stderr}"
exit 0
