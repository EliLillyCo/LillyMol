#!/bin/sh
# Wrapper for generate_profile

dir=$(dirname $0)

script="${0%%.sh}.py"
echo $script

python ${script} --feature_descriptions ${dir}/column_descriptions.txt "$@"
