#!/usr/bin/env bash

# Write to stdout a directory suitable for bazel's --output_user_root=.
# Diagnostic messages must go to stderr because callers capture stdout and pass
# it directly to bazel.

declare -i inside_lilly
if [[ $(hostname -d 2>/dev/null) =~ 'lilly.com' ]] ; then
    inside_lilly=1
else
    inside_lilly=0
fi

# bazel will not work on an NFS mounted file system. So if you are on an NFS
# file system, you must specify a value for --output_user_root that is locally
# mounted. Note that the bazel cache can get quite large, 1-2GB.

# If inside Lilly, use local scratch storage when available.
if [[ ${inside_lilly} -eq 1 && -d '/node/scratch' ]] ; then
    echo "--output_user_root=/node/scratch/${USER}/bazel"
elif [[ -d '/Projects' && -w '/Projects' ]] ; then
    echo "--output_user_root=/Projects/${USER}/bazel"
elif [[ $(uname) == 'Linux' && $(df -TP "${HOME}") =~ 'nfs' ]] ; then
    echo 'Your HOME dir is an NFS mounted file system. bazel will not work there.' >&2
    echo 'Will attempt to use /tmp for bazel cache; change this if needed.' >&2
    echo "--output_user_root=/tmp/bazel_${USER}"
else
    # If HOME is local storage, SSD recommended, bazel's default is OK.
    echo ""
fi
