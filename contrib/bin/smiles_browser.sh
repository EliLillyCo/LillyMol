#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/smiles_browser.XXXXXX")"
html="${tmpdir}/smiles_browser.html"

python3 "${SCRIPT_DIR}/smiles_browser.py" "$@" -o "${html}"

if [[ ! -s "${html}" ]] ; then
  echo "Did not create html" >&2
  exit 1
fi

if [[ "$(uname -s)" == "Darwin" ]] ; then
  open "${html}"
else
  xdg-open "${html}"
fi
