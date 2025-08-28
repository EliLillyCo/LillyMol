#!/bin/bash
here=$(dirname $0)
export PATH=$(realpath $here):$PATH
exec ruby ${here}/make_descriptors.rb "$@"
