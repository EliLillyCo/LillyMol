#!/usr/bin/env bash

# Run Python with PYTHONPATH and LD_LIBRARY_PATH pointing to the shared
# libraries staged by src/copy_shared_libraries.sh.

PYTHON=${PYTHON:-python3}

my_dir=$(dirname "$0")
soabi_file="${my_dir}/lib/lillymol.so.soabi"
actual_soabi=$(${PYTHON} -c "import sysconfig; print(sysconfig.get_config_var('SOABI') or '')")
actual_version=$(${PYTHON} -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")

if [[ -s "${soabi_file}" ]] ; then
  expected_soabi=$(tr -d '[:space:]' < "${soabi_file}")
  if [[ "${actual_soabi}" != "${expected_soabi}" ]] ; then
    echo "Python ABI mismatch for LillyMol python bindings" >&2
    echo "  ${soabi_file}: ${expected_soabi}" >&2
    echo "  ${PYTHON}: ${actual_soabi} (Python ${actual_version})" >&2
    echo "Rebuild or restage lillymol.so with the Python interpreter used by run_python.sh." >&2
    exit 1
  fi
else
  echo "Warning: cannot find ${soabi_file}; unable to verify LillyMol Python ABI" >&2
fi

if [[ -v LILLYMOL_HOME ]] ; then
  BDB="${LILLYMOL_HOME}/third_party/BDB/lib"
else
  BDB="${my_dir}/third_party/BDB/lib"
fi

if [[ ! -d ${BDB} ]] ; then
  BDB=""
fi

export LD_LIBRARY_PATH=${my_dir}/lib:${BDB}:${LD_LIBRARY_PATH}
if [[ $(uname) == 'Darwin' ]] ; then
  export DYLD_LIBRARY_PATH=${my_dir}/lib:${BDB}:${DYLD_LIBRARY_PATH}
fi
python_flags=()
if ${PYTHON} -h 2>&1 | grep -q -- " -P "; then
  python_flags=(-P)
fi

PYTHONPATH=${my_dir}/lib:${my_dir}/src:${PYTHONPATH} ${PYTHON} "${python_flags[@]}" "$@"
