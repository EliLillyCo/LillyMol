#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -r "${reactions}/butylamine_to_pyrrolidine.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -r "${reactions}/pyrrolidine_to_butylamine.rxn" -m each -o usmi -S roundtrip forward.smi
