#!/usr/bin/env bash

here=$(dirname $0)
me=$(basename $0)

exec ruby ${here}/svmfp/${me%%.sh}.rb "$@"
