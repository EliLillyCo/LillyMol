#!/usr/bin/env python

import os, sys
import logging

def most_recent_version(dir: str) ->str:
  """Scan `dir` for directories and return what looks like the most
     recent version.
     Uses python string comparisons to determine the most recent.
  """
  most_recent = ""
  for entry in os.listdir(dir):
    fname = os.path.join(dir, entry)
    if not os.path.isdir(fname):
      continue
    if entry in ['.', '..']:
      continue
    if entry > most_recent:
      most_recent = entry

  # should log error if no directories found...

  return most_recent

def get_openmp(cellar, version)->str:
  """Return the new_local_repository directive for openmp version `version`
  """
  return f"""new_local_repository(
    name = "homebrew_openmp",
    path = "{cellar}/libomp/{version}",
    build_file_content = \"\"\"
cc_library(
    name = "homebrew_openmp",
    srcs = [
       "lib/libomp.a"
    ],
    hdrs = glob(["include/*"]),
    strip_include_prefix = "include",
    visibility = ["//visibility:public"],
)
\"\"\",
)
"""

def get_gnutls(cellar, version):
  """Return the new_local_repository for gnutls version `version`.
  """
  return f"""new_local_repository(
    name = "homebrew_gnutls",
    path = "{cellar}/gnutls/{version}/",
    build_file_content = \"\"\"
cc_library(
    name = "homebrew_gnutls",
    hdrs = glob(["include/**/*.h"]),
    strip_include_prefix = "include",
    srcs = [
      "lib/libgnutls.dylib"
    ],
    visibility = ["//visibility:public"],
)
\"\"\",
)
"""

def add_homebrew_to_module_bazel(unused_args):
  """Writes homebrew repository information to stdout.
    Designed to append to an existing MODULE.bazel.
  """
  cellar = os.path.join('/opt', 'homebrew', 'Cellar')
  if not os.path.isdir(cellar):
    logging.fatal("Where is %s", cellar)

  openmp_dir = os.path.join(cellar, 'libomp')
  gnutls_dir = os.path.join(cellar, 'gnutls')

  if not os.path.isdir(openmp_dir):
    logging.fatal('Where is %s - need to install openmp', openmp_dir)

  if not os.path.isdir(gnutls_dir):
    logging.fatal('Where is %s - need to install gnutls', gnutls_dir)

  mrv = most_recent_version(openmp_dir)
  print (get_openmp(cellar, mrv))
  mrv = most_recent_version(gnutls_dir)
  print (get_gnutls(cellar, mrv))

if __name__ == '__main__':
  add_homebrew_to_module_bazel(sys.argv)
