#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -r "${reactions}/benzylamine_to_isoindoline.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -r "${reactions}/isoindoline_to_benzylamine.rxn" -m each -o usmi -S roundtrip forward.smi
