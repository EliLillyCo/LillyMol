#!/bin/bash
# Interface to random_molecular_transformations that enables
# graph edit changes only. See also minor_changes for a systematic
# approach to this problem.
# TODO: random_molecular_transformations emits too many warning
# messages, code needs cleanup..

lillymol_home=$(dirname $(dirname $0))

lib="${lillymol_home}/data/random_molecular_transformations.d"
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

# What makes this version special is the set of transformations that
# are allowed. The default file is 
probabilities="${lib}/graph_edit_changes"

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

exec random_molecular_permutations -L single=$single -L double=$double -L arom=$aromatic -L aromR=$aromatic_ring -L fsarom=$fuse_aromatic_ring -L aliphR=${aliphatic_ring} -L singleR=$single_any -L singleR=${single2} -L singleR=${single3} -r $r -R $R -P $probabilities "$@"
