#!/bin/bash

# Rename hidden files to case differentiated forms

for file in $(ls .*arom*.smi) ; do
  newname=$(echo $file | sed -e 's/^\.//' -e 's/notarom/A/g' -e 's/isarom/a/g')
  mv $file $newname
done
