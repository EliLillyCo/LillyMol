#!/bin/bash

# Rename hidden files to case insensitive forms

for file in $(ls .*arom*.smi) ; do
  newname=$(echo $file | sed -e 's/^\.//' -e 's/notarom/Al/g' -e 's/isarom/Ar/g')
  mv $file $newname
done
