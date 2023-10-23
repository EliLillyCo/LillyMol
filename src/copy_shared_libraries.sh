#!/bin/bash

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
bazel-bin/pybind/lillymol_query.so
bazel-bin/pybind/lillymol_reaction.so
bazel-bin/pybind/lillymol_standardise.so
bazel-bin/pybind/lillymol_tools.so
bazel-bin/pybind/lillymol_tsubstructure.so
bazel-bin/Molecule_Lib/*.so
)

for lib in "${libs[@]}" ; do
  if compgen -G "${lib}" > /dev/null; then
    echo "Copying ${lib}"
    cp -f ${lib} ${destdir}
  else
    echo "${lib} not found"
  fi
done

# Copy python generated protos

declare -a pb2=(
  atom_type_ext_pb2.py
  geometric_constraints_pb2.py
  mol2graph_pb2.py
  substructure_pb2.py
  reaction_pb2.py
  toggle_kekule_form_pb2.py
)

for proto_pb2 in "${pb2[@]}" ; do
  echo "Copying ${proto_pb2}"
  cp -f bazel-bin/Molecule_Lib/${proto_pb2} Molecule_Lib
done
