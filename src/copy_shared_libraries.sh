#!/usr/bin/env bash

# when building a pybind environment, we need to copy the shared libraries generated
# by bazel to somewhere permanent.
# we also copy _pb2.py files

if [[ -z $1 ]] ; then
  destdir='/lrlhps/users/rx87690/LillyMolPrivate/lib'
else
  destdir=$1
fi

if [[ ! -d ${destdir} ]] ; then
  mkdir -p ${destdir}
fi

declare -a libs=(
bazel-bin/pybind/lillymol.so
bazel-bin/pybind/lillymol_io.so
bazel-bin/pybind/lillymol_fingerprint.so
bazel-bin/pybind/lillymol_query.so
bazel-bin/pybind/lillymol_reaction.so
bazel-bin/pybind/lillymol_standardise.so
bazel-bin/pybind/lillymol_tools.so
bazel-bin/pybind/lillymol_tsubstructure.so
bazel-bin/Molecule_Lib/*.so
)

for lib in "${libs[@]}" ; do
  if compgen -G "${lib}" > /dev/null; then
    libname=$(basename ${lib})
    if [[ ! -s ${destdir}/${libname} || ${lib} -nt ${destdir}/${libname} ]] ; then
      echo "Copying ${lib}"
      cp -f ${lib} ${destdir}
    fi
  else
    echo "${lib} not found"
  fi
done

# Copy python generated protos

declare -a pb2=(
  atom_type_ext_pb2.py
  geometric_constraints_pb2.py
  mol2graph_pb2.py
  pharmacophore_pb2.py
  substructure_pb2.py
  reaction_pb2.py
  toggle_kekule_form_pb2.py
)

for proto_pb2 in "${pb2[@]}" ; do
  source="bazel-bin/Molecule_Lib/${proto_pb2}"
  dest="Molecule_Lib/${proto_pb2}"

  if [[ ! -s ${dest} || ${source} -nt ${dest} ]] ; then
    echo "Copying ${proto_pb2}"
    cp -f bazel-bin/Molecule_Lib/${proto_pb2} Molecule_Lib
  fi
done

declare -a tools_pb2=(
  dicer_fragments_pb2.py
  iwdescr_pb2.py
  xlogp_pb2.py
)
for proto_pb2 in "${tools_pb2[@]}" ; do
  source="bazel-bin/Molecule_Tools/${proto_pb2}"
  dest="Molecule_Tools/${proto_pb2}"
  if [[ ! -s ${dest} || ${source} -nt ${dest} ]] ; then
    echo "Copying ${proto_pb2}"
    cp -f bazel-bin/Molecule_Tools/${proto_pb2} Molecule_Tools
  fi
done
