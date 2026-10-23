#!/usr/bin/env bash
here=$(dirname $0)
exec ruby ${here}/molecule_filter_parallel.rb "$@"
