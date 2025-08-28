# Used during Docker builds to update WORKSPACE to
# reflect the python installation in the Docker file
# Basically just finds the relevant 'path = ' sections
# in WORKSPACE and updates them
# We assume that the caller has copied the original to a
# safe location, this script writes to stdout.
#    cp WORKSPACE /tmp
#    update_python_in_workspace /tmp/WORKSPACE > WORKSPACE

import logging
import os
import re
import sys
import sysconfig

have_pybind = True
try:
  from pybind11 import *
except:
  logging.info("No pybind11 build will fail")
  have_pybind = False

# Scan a WORKSPACE file and change the 'path = ' directive
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
  update_python_in_workspace(sys.argv)
