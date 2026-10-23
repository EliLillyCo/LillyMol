#!/bin/bash

set -e

nplotnn=$1
datadir=$4

bindir=$(dirname "$nplotnn")
nn2proto="$bindir/nn2proto"
expected=$(dirname "$0")/../default/out/stdout

if [ ! -x "$nn2proto" ]
then
  echo "Cannot find nn2proto in $bindir" >&2
  exit 1
fi

$nn2proto -T nn1.tfrecord "$datadir/nn1.nn"
$nplotnn -X tfdata nn1.tfrecord > tfdata.out
diff -w tfdata.out "$expected"
echo "TFDataRecord nplotnn output matches text output"
