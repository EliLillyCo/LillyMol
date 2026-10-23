#!/bin/bash

# Rename hidden files to case insensitive forms

# If it looks like the files are already in the correct form, do nothing

if [[ -s 'rings_6Al.smi' ]] ; then
  exit 0
fi

for file in $(ls rings_*.smi) ; do
  newname=$(echo $file | sed -e 's/A/Al/g' -e 's/a/Ar/g')
  mv $file $newname
done
