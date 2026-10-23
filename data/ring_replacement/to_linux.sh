#!/bin/bash

# Rename files from Mac forms to Linux forms.

# If it looks like the files are already in the correct form, do nothing.
if [[ -s rings_6a.smi ]] ; then
  exit 0
fi

for file in $(ls rings*.smi) ; do
  newname=$(echo $file | sed -e 's/Ar/a/g' -e 's/Al/A/g')
  mv $file $newname
done
