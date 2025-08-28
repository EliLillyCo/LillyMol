#!/bin/bash
here=$(dirname $(realpath $0))

me=$(basename $0)
exec ruby ${here}/${me/.sh/.rb} "$@"
