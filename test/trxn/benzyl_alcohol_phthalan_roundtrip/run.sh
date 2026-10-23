#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -r "${reactions}/benzyl_alcohol_to_phthalan.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -r "${reactions}/phthalan_to_benzyl_alcohol.rxn" -m each -o usmi -S roundtrip forward.smi
