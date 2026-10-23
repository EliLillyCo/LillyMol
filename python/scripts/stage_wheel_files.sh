#!/usr/bin/env bash
set -euo pipefail

# Stage LillyMol nanobind extension modules and runtime Python files for wheel
# building. This script assumes copy_shared_libraries.sh has already populated
# ${LILLYMOL_HOME}/lib. It does not alter ${LILLYMOL_HOME}/lib or run_python.sh.

script_dir=$(cd "$(dirname "$0")" && pwd)
python_dir=$(cd "${script_dir}/.." && pwd)
repo_home=$(cd "${python_dir}/.." && pwd)
src_dir="${repo_home}/src"
lib_dir="${LILLYMOL_LIB_DIR:-${repo_home}/lib}"
stage_dir="${1:-${python_dir}/prebuilt}"

if [[ ! -d "${lib_dir}" ]] ; then
  echo "Library directory not found: ${lib_dir}" >&2
  echo "Run src/copy_shared_libraries.sh ${lib_dir} first." >&2
  exit 1
fi

rm -rf "${stage_dir}"
mkdir -p "${stage_dir}"

nanobind_modules=(
  lillymol
  lillymol_bdb
  lillymol_gfp_server
)

for module in "${nanobind_modules[@]}" ; do
  source="${lib_dir}/${module}.so"
  if [[ ! -s "${source}" ]] ; then
    fallback="${src_dir}/bazel-bin/nanobind/${module}.so"
    if [[ -s "${fallback}" ]] ; then
      source="${fallback}"
    else
      echo "Missing nanobind module ${module}.so in ${lib_dir} and bazel-bin/nanobind" >&2
      exit 1
    fi
  fi
  cp -f "${source}" "${stage_dir}/"
done

# Compatibility modules for older split-module imports are pure Python shims in
# the nanobind build. They re-export symbols from the canonical lillymol module.
for shim in "${src_dir}"/nanobind/python_shims/*.py ; do
  if [[ -s "${shim}" ]] ; then
    cp -f "${shim}" "${stage_dir}/"
  fi
done

# Private shared libraries needed by extension modules. Keep these at wheel top
# level; stage copies can later be patched with an $ORIGIN rpath if needed.
while IFS= read -r lib ; do
  cp -f "${lib}" "${stage_dir}/"
done < <(find "${lib_dir}" -maxdepth 1 -type f -name 'lib*.so*' | sort)

copy_python_file() {
  local source="$1"
  local destination="$2"
  if [[ -s "${source}" ]] ; then
    mkdir -p "$(dirname "${destination}")"
    cp -f "${source}" "${destination}"
  fi
}

copy_generated_python_file() {
  local relative_path="$1"
  local destination="${stage_dir}/${relative_path}"

  if [[ -s "${src_dir}/${relative_path}" ]] ; then
    copy_python_file "${src_dir}/${relative_path}" "${destination}"
  elif [[ -s "${src_dir}/bazel-bin/${relative_path}" ]] ; then
    copy_python_file "${src_dir}/bazel-bin/${relative_path}" "${destination}"
  fi
}

# Python protobuf modules used by existing bindings/tests and GFP HTTP helpers.
for proto in \
  atom_type_ext_pb2.py \
  geometric_constraints_pb2.py \
  mol2graph_pb2.py \
  molecule_to_query_pb2.py \
  pharmacophore_pb2.py \
  substructure_pb2.py \
  reaction_pb2.py \
  toggle_kekule_form_pb2.py ; do
  copy_generated_python_file "Molecule_Lib/${proto}"
done

for proto in \
  common_names_pb2.py \
  dicer_fragments_pb2.py \
  iwdescr_pb2.py \
  xlogp_pb2.py ; do
  copy_generated_python_file "Molecule_Tools/${proto}"
done

for proto in \
  nearneighbours_pb2.py \
  nn_request_pb2.py ; do
  copy_generated_python_file "Utilities/GFP_Tools/${proto}"
done

copy_python_file "${src_dir}/Utilities/GFP_Tools/gfp_http_server.py" "${stage_dir}/Utilities/GFP_Tools/gfp_http_server.py"
copy_python_file "${src_dir}/Utilities/GFP_Tools/gfp_client_app.py" "${stage_dir}/Utilities/GFP_Tools/gfp_client_app.py"

# Make package directories explicit. Namespace packages would work, but explicit
# __init__.py files make the wheel contents easier to inspect and debug.
for package_dir in \
  "${stage_dir}/Molecule_Lib" \
  "${stage_dir}/Molecule_Tools" \
  "${stage_dir}/Utilities" \
  "${stage_dir}/Utilities/GFP_Tools" ; do
  if [[ -d "${package_dir}" ]] ; then
    touch "${package_dir}/__init__.py"
  fi
done

# If available, make staged shared objects look for private LillyMol shared
# libraries beside themselves after wheel installation. This modifies only the
# staged copies, not ${LILLYMOL_HOME}/lib.
if command -v patchelf >/dev/null 2>&1 ; then
  while IFS= read -r sofile ; do
    chmod u+w "${sofile}" || true
    patchelf --set-rpath '$ORIGIN' "${sofile}"
  done < <(find "${stage_dir}" -maxdepth 1 -type f -name '*.so*' | sort)
else
  echo "patchelf not found; wheel may require LD_LIBRARY_PATH or auditwheel repair" >&2
fi

cat <<EOF
Staged LillyMol nanobind wheel files in ${stage_dir}

Next steps:
  cd ${python_dir}
  python -m pip install build wheel setuptools
  python -m build --wheel

  or

  uv pip install build wheel setuptools
  uv build --wheel
EOF
