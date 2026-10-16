#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -P "${reactions}/phenethylamine_to_indoline.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -P "${reactions}/indoline_to_phenethylamine.rxn" -m each -o usmi -S roundtrip forward.smi
