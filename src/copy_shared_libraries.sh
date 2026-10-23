#!/usr/bin/env bash

# Copy nanobind extension modules, private LillyMol shared-library dependencies,
# and generated Python protos to a directory suitable for PYTHONPATH.
#
# Nanobind is the production Python binding. The older pybind tree remains in
# the repository for now, but this script no longer stages pybind extension
# modules.

set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo_home=$(cd "${here}/.." && pwd)

if [[ $# -gt 0 ]] ; then
  destdir=$1
else
  destdir="${repo_home}/lib"
fi

mkdir -p "${destdir}"

# Remove stale pybind-era extension modules that used to be staged in lib. A
# stale lillymol.so would shadow the newly copied nanobind lillymol.so; stale
# split modules are worse because they can import old pybind11 types into a
# nanobind process.
rm -f \
  "${destdir}/lillymol_atom.so" \
  "${destdir}/lillymol_bond.so" \
  "${destdir}/lillymol_fingerprint.so" \
  "${destdir}/lillymol_io.so" \
  "${destdir}/lillymol_query.so" \
  "${destdir}/lillymol_reaction.so" \
  "${destdir}/lillymol_ring.so" \
  "${destdir}/lillymol_set_of_atoms.so" \
  "${destdir}/lillymol_standardise.so" \
  "${destdir}/lillymol_tools.so" \
  "${destdir}/lillymol_tsubstructure.so" \
  "${destdir}/lillymol.so" \
  "${destdir}/lillymol.abi3.so" \
  "${destdir}/lillymol_bdb.so" \
  "${destdir}/lillymol_bdb.abi3.so" \
  "${destdir}/lillymol_depict.so" \
  "${destdir}/lillymol_depict.abi3.so" \
  "${destdir}/lillymol_gfp_server.so" \
  "${destdir}/lillymol_gfp_server.abi3.so" \
  "${destdir}/lillymol_nb.so" \
  "${destdir}/lillymol_nb.abi3.so" \
  "${destdir}/lillymol_nb_bdb.so" \
  "${destdir}/lillymol_nb_bdb.abi3.so" \
  "${destdir}/lillymol.so.soabi"

copy_if_newer() {
  local source=$1
  local libname
  libname=$(basename "${source}")

  if [[ ! -e "${source}" ]] ; then
    return
  fi

  if [[ ! -s "${destdir}/${libname}" || "${source}" -nt "${destdir}/${libname}" ]] ; then
    echo "Copying ${source}"
    cp -f "${source}" "${destdir}"
  fi
}

copy_required_glob() {
  local pattern=$1
  local found=0

  shopt -s nullglob
  for lib in ${pattern}; do
    found=1
    copy_if_newer "${lib}"
  done
  shopt -u nullglob

  if [[ ${found} -eq 0 ]] ; then
    echo "${pattern} not found" >&2
    return 1
  fi
}

copy_optional_glob() {
  local pattern=$1

  shopt -s nullglob
  for lib in ${pattern}; do
    copy_if_newer "${lib}"
  done
  shopt -u nullglob
}

copy_python_proto_if_newer() {
  local source=$1
  local dest=$2

  if [[ ! -s "${source}" ]] ; then
    return
  fi
  if [[ ! -s "${dest}" || "${source}" -nt "${dest}" ]] ; then
    echo "Copying $(basename "${source}")"
    mkdir -p "$(dirname "${dest}")"
    cp -f "${source}" "${dest}"
  fi
}

cd "${here}"

copy_if_newer "bazel-bin/nanobind/lillymol.so"
copy_optional_glob "bazel-bin/nanobind/lillymol_bdb.so"
copy_optional_glob "bazel-bin/nanobind/lillymol_depict.so"
copy_optional_glob "bazel-bin/nanobind/lillymol_gfp_server.so"

copy_if_newer "bazel-bin/nanobind/lillymol.so.soabi"

copy_optional_glob "bazel-bin/Molecule_Lib/lib*.so"
copy_optional_glob "bazel-bin/Molecule_Tools/lib*.so"
copy_optional_glob "bazel-bin/Molecule_Tools_Bdb/lib*.so"
copy_optional_glob "bazel-bin/Utilities/GFP_Tools/lib*.so"

for shim in nanobind/python_shims/*.py; do
  copy_if_newer "${shim}"
done

for proto_pb2 in \
  atom_type_ext_pb2.py \
  geometric_constraints_pb2.py \
  mol2graph_pb2.py \
  molecule_to_query_pb2.py \
  pharmacophore_pb2.py \
  standardise_pb2.py \
  substructure_pb2.py \
  reaction_pb2.py \
  toggle_kekule_form_pb2.py; do
  copy_python_proto_if_newer "bazel-bin/Molecule_Lib/${proto_pb2}" "Molecule_Lib/${proto_pb2}"
done

for proto_pb2 in \
  common_names_pb2.py \
  dicer_fragments_pb2.py \
  iwdescr_pb2.py \
  xlogp_pb2.py; do
  copy_python_proto_if_newer "bazel-bin/Molecule_Tools/${proto_pb2}" "Molecule_Tools/${proto_pb2}"
done

for proto_pb2 in \
  nearneighbours_pb2.py \
  nn_request_pb2.py; do
  copy_python_proto_if_newer "bazel-bin/Utilities/GFP_Tools/${proto_pb2}" "Utilities/GFP_Tools/${proto_pb2}"
done
