#!/usr/bin/env bash

# Parts of LillyMol come in compressed form or with obfuscated file names.
# Convert those.
# If this script has been run and succeeded it will be a no-op

dicer_fragments_dir="${LILLYMOL_HOME}/data/dicer_fragments/"
if [[ -s "${dicer_fragments_dir}/dicer_fragments.tar.xz" ]] ; then
  echo "Uncompressing fragments"
  (cd ${dicer_fragments_dir} && tar Jxf dicer_fragments.tar.xz)
fi

chembl_sidechains="${LILLYMOL_HOME}/data/chembl_sidechains.textproto.xz"
if [[ -s ${chembl_sidechains} ]] ; then
  echo "Uncompressing CHEMBL sidechains"
  xz --decompress ${chembl_sidechains}
fi

ring_replacement="${LILLYMOL_HOME}/data/ring_replacement/"

if [[ -d "${ring_replacement}" ]] ; then
  if [[ $(uname) == "Linux" ]] ; then
    (cd ${LILLYMOL_HOME}/data/ring_replacement && ./to_linux.sh)
  fi
fi
