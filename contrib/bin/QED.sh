#!/bin/bash
here=$(dirname $0)

if [[ -v LILLYMOL_HOME ]] ; then
  qed_dir="${LILLYMOL_HOME}/data/QED"
else
  qed_dir="${here}/../../data/QED"
fi

if [[ ! -d ${qed_dir} ]] ; then
  echo "Where is my QED data '${qed_dir}'" >&2
  exit 1
fi

./bazel-bin/Molecule_Tools/qed -Q QUERY=S:${qed_dir}/unwanted-groups.smt \
        -Q ACC=S:${qed_dir}/acceptor.smt \
        -Q DON=S:${qed_dir}/donor.smt \
        -Q ROTB=S:${qed_dir}/rotb.smt "$@"
