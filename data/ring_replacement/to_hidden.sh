#!/bin/bash
# Case insensitive file systems cannot differentiate 5a5A.smi from 5a5a.smi.
# Convert all ring replacement files in the current directory to case differentiated forms.

for file in $(/bin/ls rings_*[3-9][a,A]*.smi) ; do 
  newname=$(echo ${file} | sed -E -e 's/([3-9])A/\1notarom/g' -e 's/([3-9])a/\1isarom/g')
  mv ${file} .${newname}
done
