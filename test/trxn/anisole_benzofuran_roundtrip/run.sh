#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -P "${reactions}/anisole_to_benzofuran.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -P "${reactions}/benzofuran_to_anisole.rxn" -m each -o usmi -S roundtrip forward.smi
