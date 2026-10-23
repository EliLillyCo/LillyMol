#!/usr/bin/env bash

# LillyMol likely needs libraries from HomeBrew.
# Append the likely location to stdin.
# add_homebrew_to_modele_bazel.sh >> MODULE.bazel 

cellar='/opt/homebrew/Cellar'
cellar='/tmp/q'

if [[ ! -d ${cellar} ]] ; then
  echo "$0: No ${cellar} cannot continue" >&2
  exit 1
fi

if [[ ! "${cellar}/libomp" ]] ; then
  echo "libomp not found in ${cellar}" >&2
  exit 1
fi

if [[ ! "${cellar}/gnutls" ]] ; then
  echo "gnutls not found in ${cellar}" >&2
  exit 1
fi

# Return the most recent version of a package in the homebrew cellar
function most_recent_version {
  dir=$1
  most_recent=""
  for fname in "${dir}/*" ; do
    if [[ ${most_recent} == "" ]] ; then
      most_recent=${fname}
    elif [[ ${fname} > ${most_recent} ]] ; then
      most_recent=${fname}
    fi
  done
  echo $most_recent
}

openmp=$(cat << EOF
new_local_repository(
    name = "homebrew_openmp",
    path = "/opt/homebrew/Cellar/libomp/20.1.5/",
    build_file_content = """
cc_library(
    name = "homebrew_openmp",
    srcs = [
       "lib/libomp.a"
    ],
    hdrs = glob(["include/*"]),
    strip_include_prefix = "include",
    visibility = ["//visibility:public"],
)
""",
)
EOF
)

gnutls=$(cat << 'EOF'
new_local_repository(
    name = "homebrew_gnutls",
    path = "/opt/homebrew/Cellar/gnutls/3.8.9/",
    build_file_content = """
cc_library(
    name = "homebrew_gnutls",
    hdrs = glob(["include/**/*.h"]),
    strip_include_prefix = "include",
    srcs = [
      "lib/libgnutls.dylib"
    ],
    visibility = ["//visibility:public"],
)
""",
)
EOF
)
g=most_recent_version("cellar/openmp")
echo $g
