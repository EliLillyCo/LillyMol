#!/usr/bin/env bash

here=$(dirname $0)
exec ruby ${here}/svmfp/svmfp_make.rb "$@"
