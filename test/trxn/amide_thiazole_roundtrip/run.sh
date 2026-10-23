#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -P "${reactions}/amide_to_thiazole.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -P "${reactions}/thiazole_to_amide.rxn" -m each -o usmi -S roundtrip forward.smi
