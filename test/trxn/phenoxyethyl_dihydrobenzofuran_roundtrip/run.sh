#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -r "${reactions}/phenoxyethyl_to_dihydrobenzofuran.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -r "${reactions}/dihydrobenzofuran_to_phenoxyethyl.rxn" -m each -o usmi -S roundtrip forward.smi
