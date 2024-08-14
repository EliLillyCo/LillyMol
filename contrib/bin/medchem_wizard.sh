#!/bin/bash

if [ -v LILLYMOL_HOME ] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(readlink -e $0))))
fi

rxn_dir=${LILLYMOL_HOME}/data/MedchemWizard

if [[ ! -d "${rxn_dir}" ]] ; then
  echo "$0 where are my reactions ${rxn_dir}" >&2
  exit 1
fi

exe=${LILLYMOL_HOME}/bin/$(uname)/medchemwizard

exec ${exe} -R "${rxn_dir}/REACTIONS" $@
