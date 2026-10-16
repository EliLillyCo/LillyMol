#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -r "${reactions}/biphenyl_to_fluorene.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -r "${reactions}/fluorene_to_biphenyl.rxn" -m each -o usmi -S roundtrip forward.smi
