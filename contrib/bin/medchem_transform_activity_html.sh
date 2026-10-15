#!/usr/bin/env bash
# Use the staged LillyMol Python ABI and libraries.
set -e
script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export LILLYMOL_HOME=${LILLYMOL_HOME:-"${script_dir}/../.."}
exec "${LILLYMOL_HOME}/run_python.sh" "${script_dir}/medchem_transform_activity_html.py" "$@"
