#!/bin/bash

if [[ -v LILLYMOL_HOME ]] ; then
  true
else
  export LILLYMOL_HOME=$(dirname $(dirname $(dirname $(readlink -e $0))))
fi

lib="${LILLYMOL_HOME}/data/random_molecular_permutations.d"
if [[ ! -d ${lib} ]] ; then
  echo "Missing data ${lib}" >&2
  ls ${lib}
  exit 1
fi

# Define transformation libraries
single="${lib}/single_attachment_fragments.smi";
single_any="${lib}/single_attachment_fragments_any_atom.smi"
double="${lib}/double_attachment_fragments.smi";
aromatic="${lib}/aromatic_attachment_fragments.smi";
aromatic_ring="${lib}/aromatic_ring.smi"
fuse_aromatic_ring="${lib}/fuse_aromatic_ring.smi"
aliphatic_ring="${lib}/aliphatic_ring.smi"

single2="${lib}/Fragments.rings_aliphatic.smi"
single3="${lib}/Fragments.rings_aromatic.smi"

probabilities="${lib}/random_molecular_transformations"

# default min and max ring size, can be overwritten
# in the command r and/or R will get repeated, but that does not matter
    r=4
    R=7
    for ((i=1; i<$#; i++)); do
        if [[ "${!i}" == "-r" ]]; then
            j=$((i+1))
            r="${!j}"
        elif [[ "${!i}" == "-R" ]]; then
            j=$((i+1))
            R="${!j}"
        fi
    done

exe=${LILLYMOL_HOME}/bin/$(uname)/random_molecular_permutations
exec ${exe} -L single=$single -L double=$double -L arom=$aromatic -L aromR=$aromatic_ring -L fsarom=$fuse_aromatic_ring -L aliphR=${aliphatic_ring} -L singleR=$single_any -L singleR=${single2} -L singleR=${single3} -r $r -R $R -P $probabilities "$@"
