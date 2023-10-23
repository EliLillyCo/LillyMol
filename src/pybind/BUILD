
package(default_visibility = ["//visibility:public"])

load("@rules_python_pytest//python_pytest:defs.bzl", "py_pytest_test")

cc_library(
    name = "lillymol_pybind",
    srcs = [
        "lillymol_pybind.cc",
    ],
    deps = [
        ":molecule",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
    copts = [
        "-save-temps",
    ],
)

cc_shared_library(
    name = "lillymol",
    deps = [
        ":lillymol_pybind",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
  ],
  shared_lib_name = "lillymol.so",
)

cc_library(
    name = "lillymol_pybind_substructure_query",
    srcs = [
        "lillymol_pybind_query.cc",
    ],
    deps = [
        ":molecule",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
    name = "lillymol_query",
    deps = [
        ":lillymol_pybind_substructure_query",
    ],
    shared_lib_name = "lillymol_query.so",
)

cc_library(
    name = "lillymol_pybind_chemical_standardisation",
    srcs = [
        "lillymol_pybind_std.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
  name = "lillymol_standardise",
  deps = [
      ":lillymol_pybind_chemical_standardisation",
  ],
  shared_lib_name = "lillymol_standardise.so",
)

cc_library(
    name = "lillymol_pybind_set_of_atoms",
    srcs = [
        "lillymol_pybind_set_of_atoms.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
  name = "lillymol_set_of_atoms",
  deps = [
      ":lillymol_pybind_set_of_atoms",
  ],
  shared_lib_name = "lillymol_set_of_atoms.so",
)

cc_library(
    name = "lillymol_pybind_ring",
    srcs = [
        "lillymol_pybind_ring.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
  name = "lillymol_ring",
  deps = [
      ":lillymol_pybind_ring",
  ],
  shared_lib_name = "lillymol_ring.so",
)

cc_library(
    name = "lillymol_pybind_atom",
    srcs = [
        "lillymol_pybind_atom.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
  name = "lillymol_atom",
  deps = [
      ":lillymol_pybind_atom",
  ],
  shared_lib_name = "lillymol_atom.so",
)

cc_library(
    name = "lillymol_pybind_bond",
    srcs = [
        "lillymol_pybind_bond.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
  name = "lillymol_bond",
  deps = [
      ":lillymol_pybind_bond",
  ],
  shared_lib_name = "lillymol_bond.so",
)

#cc_library(
#    name = "lillymol_pybind_chiral_centre",
#    srcs = [
#        "lillymol_pybind_chiral_centre.cc",
#    ],
#    deps = [
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:moleculeio",
#        "@pybind11",
#        "@python",
#    ],
#)

#cc_shared_library(
#  name = "lillymol_chiral_centre",
#  deps = [
#      ":lillymol_pybind_chiral_centre",
#  ],
#  shared_lib_name = "lillymol_chiral_centre.so",
#)

cc_library(
    name = "lillymol_pybind_io",
    srcs = [
        "lillymol_pybind_io.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
  name = "lillymol_io",
  deps = [
      ":lillymol_pybind_io",
  ],
  shared_lib_name = "lillymol_io.so",
)

cc_library(
    name = "molecule",
    srcs = [
        "molecule.cc",
    ],
    hdrs = [
        "molecule.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_shared_library(
  name = "molecule_shared",
  deps = [
      ":molecule",
  ],
  shared_lib_name = "lillymol_molecule_misc.so",
)

cc_library(
    name = "lillymol_pybind_reaction",
    srcs = [
        "lillymol_pybind_reaction.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "@pybind11",
        "@python",
    ],
)

cc_shared_library(
  name = "lillymol_reaction",
  deps = [
      ":lillymol_pybind_reaction",
  ],
  shared_lib_name = "lillymol_reaction.so",
)

cc_library(
    name = "tsubstructure",
    srcs = [
        "tsubstructure.cc",
	"lillymol_pybind_tsubstructure.cc",
    ],
    hdrs = [
        "tsubstructure.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
	"@pybind11",
        "@python",
  ],
)
cc_shared_library( 
  name = "tsubstructure_shared",
  deps = [
      ":tsubstructure",
  ],
  shared_lib_name = "lillymol_tsubstructure.so",
)

cc_library(
    name = "lillymol_tools",
    srcs = [
        "lillymol_pybind_tools.cc",
    ],
    hdrs = [
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:unique_molecules_api",
	"@pybind11",
        "@python",
    ],
)

cc_shared_library(
    name = "lillymol_tools_shared",
    deps = [
        ":lillymol_tools",
    ],
    shared_lib_name = "lillymol_tools.so",
)

# Cannot figure out how to make this work. TODO:ianwatson
#py_pytest_test(
#    name = "test_lillymol",
#    size = "small",
#    srcs = ["lillymol_test.py"],
#    deps = [
#      ":lillymol_pybind",
#    ],
#)

