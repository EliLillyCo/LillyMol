#!/bin/bash
set -euo pipefail
trxn=$1
indir=$2
reactions="${LILLYMOL_HOME}/data/MedchemWizard"
"${trxn}" -r "${reactions}/phenylpropylamine_to_tetrahydroquinoline.rxn" -m each -o usmi -S forward "${indir}/input.smi"
"${trxn}" -r "${reactions}/tetrahydroquinoline_to_phenylpropylamine.rxn" -m each -o usmi -S roundtrip forward.smi
