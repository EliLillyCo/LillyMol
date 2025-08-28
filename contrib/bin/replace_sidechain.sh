#!/bin/bash

here=$(dirname $0)

exec ruby ${here}/replace_sidechain.rb "$@"
