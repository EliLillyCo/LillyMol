#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -P "${reactions}/amide_to_1_2_4_oxadiazole.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -P "${reactions}/1_2_4_oxadiazole_to_amide.rxn" -m each -o usmi -S roundtrip forward.smi
