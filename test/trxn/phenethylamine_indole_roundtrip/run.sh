#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -P "${reactions}/phenethylamine_to_indole.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -P "${reactions}/indole_to_phenethylamine.rxn" -m each -o usmi -S roundtrip forward.smi
