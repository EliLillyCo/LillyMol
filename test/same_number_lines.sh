#!/usr/bin/env bash

# Returns 0 if the two file arguments have the same line count.

if [[ $# -ne 2 ]] ; then
  echo "Must specify two files" >&2
  exit 1
fi

l1=$(wc -l < $1)
l2=$(wc -l < $2)

if [[ ${l1} -eq ${l2} ]] ; then
  exit 0
else
  exit 1
fi
