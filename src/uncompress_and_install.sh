#!/bin/bash

# Parts of LillyMol come in compressed form or with obfuscated file names.
# Convert those.
# If thi script has been run and succeeded it will be a no-op

dicer_fragments_dir="${LILLYMOL_HOME}/data/dicer_fragments/"
if [[ -s "${dicer_fragments_dir}/dicer_fragments.tar.xz" ]] ; then
  echo "Uncompressing fragments"
  cd ${dicer_fragments_dir} && tar Jxf dicer_fragments.tar.xz
fi

hidden_ring="${LILLYMOL_HOME}/data/ring_replacement/.rings_6isarom.smi"

if [[ -s "${hidden_ring}" ]] ; then
  if [[ $(uname) == "Linux" ]] ; then
    echo "Unhiding replacement rings"
    cd ${LILLYMOL_HOME}/data/ring_replacement && ./to_linux.sh
  fi
  # not doing mac yet
fi
