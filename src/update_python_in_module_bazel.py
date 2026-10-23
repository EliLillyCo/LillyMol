# Used during Docker builds to update MODULE.bazel to
# reflect the python installation in the Docker file
# Basically just finds the relevant 'path = ' sections
# in MODULE.bazel and updates them
# We assume that the caller has copied the original to a
# safe location, this script writes to stdout.
#    cp MODULE.bazel /tmp
#    update_python_in_workspace /tmp/MODULE.bazel > MODULE.bazel

import os
import re
import sys
import sysconfig

from absl import app
from absl import logging

have_pybind = True
try:
  from pybind11 import *
except:
  logging.info("No pybind11 build will fail")
  have_pybind = False

# From Gemini, much better than what I had.

def update_bazel_python_path(argv):
    # The location of python and pybind11 includes
    python_install = sysconfig.get_path('include')
    logging.info("Python in %s", python_install)

    pybind_install = ""
    if have_pybind:
        pybind_install = get_include()
        logging.info("pybind11 in %s", pybind_install)

    with open(argv[1], 'r') as f:
        content = f.read()

    # Regex breakdown:
    # 1. Finds 'new_local_repository'
    # 2. Ensures 'name = "python"' is the next argument
    # 3. Matches 'path = "..."' and captures the groups around the actual path string
    pattern = r'(new_local_repository\(\s*name\s*=\s*"python",.*?path\s*=\s*")([^"]*)(")'
    
    # re.DOTALL is critical here so the '.' matches newlines between 'name' and 'path'
    new_content = re.sub(pattern, r'\1' + python_install + r'\3', content, flags=re.DOTALL)

    if new_content == content:
        logging.error("No changes made. Check if the 'name' and 'path' structure matches.")
    else:
        print(new_content)

# Scan a MODULE.bazel file and change the 'path = ' directive
# for new_local_repository's 'pybind11' and 'python'.
# We assume that each new_local_repository contains
# 'name = ' BEFORE 'path ='.

def update_python_in_workspace(argv):
  # The location of python and pybind11 includes
  python_install = sysconfig.get_path('include')
  logging.info("Python in %s", python_install)

  pybind_install = ""
  if have_pybind:
    pybind_install = get_include()
    logging.info("pybind11 in %s", pybind_install)

  in_new_local_repository = False
  # The name of the new_local_repository we are in
  name = ""

  with open(argv[1], "r") as reader:
    for line in reader:
      line = line.rstrip()
      if line.startswith("new_local_repository"):
        in_new_local_repository = True
        print(line)
        continue

      if line.startswith(")"):
        if in_new_local_repository:
          name = ""
          in_new_local_repository = False
        print(line)
        continue

      m = re.search(r'name *= *"(\S+)"', line)
      if m:
        name = m[1]
        print(line)
        continue
      
      m = re.search(r'path *= *"(\S+)"', line)
      if m:
        if not in_new_local_repository:
          print(line)
        elif name == "python":
          print('    path = "' + python_install + '",')
        elif name == "pybind11" and len(pybind_install) > 0:
          print('    path = "' + pybind_install + '",')
        else:
          print(line)
      else:
        print(line)

if __name__ == "__main__":
  app.run(update_bazel_python_path)
