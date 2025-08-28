# load("@rules_python_pytest//python_pytest:defs.bzl", "py_pytest_test")
# load("@rules_python//python:py_test.bzl", "py_test")
load("@pybind11_bazel//:build_defs.bzl", "pybind_extension")
load("@rules_python//python:defs.bzl", "py_library", "py_test")

package(default_visibility = ["//visibility:public"])

pybind_extension(
    name = "lillymol",
    srcs = [
        "lillymol_pybind.cc",
    ],
    copts = [
        "-save-temps",
    ],
    deps = [
        ":molecule",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:alogp_lib",
        "//Molecule_Tools:xlogp_lib",
        "@python",
#       "@third_party//:python",
    ],
)

pybind_extension(
    name = "lillymol_fingerprint",
    srcs = [
        "lillymol_pybind_fingerprint.cc",
    ],
    deps = [
        ":molecule",
        "//Molecule_Lib:atom_pair_fingerprint",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:linear_fingerprint",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:ec_fingerprint_lib",
        "@python",
    ],
)


pybind_extension(
    name = "lillymol_query",
    srcs = [
        "lillymol_pybind_query.cc",
    ],
    deps = [
        ":molecule",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
)

pybind_extension(
    name = "lillymol_standardise",
    srcs = [
        "lillymol_pybind_std.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
)

pybind_extension(
    name = "lillymol_set_of_atoms",
    srcs = [
        "lillymol_pybind_set_of_atoms.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
)

pybind_extension(
    name = "lillymol_ring",
    srcs = [
        "lillymol_pybind_ring.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
)

pybind_extension(
    name = "lillymol_atom",
    srcs = [
        "lillymol_pybind_atom.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
)

pybind_extension(
    name = "lillymol_bond",
    srcs = [
        "lillymol_pybind_bond.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
)

pybind_extension(
    name = "lillymol_io",
    srcs = [
        "lillymol_pybind_io.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
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

pybind_extension(
    name = "lillymol_reaction",
    srcs = [
        "lillymol_pybind_reaction.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "@python",
    ],
)

pybind_extension(
    name = "lillymol_tsubstructure",
    srcs = [
        "lillymol_pybind_tsubstructure.cc",
    ],
    deps = [
        ":tsubstructure",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_library(
    name = "tsubstructure",
    srcs = [
        "tsubstructure.cc",
    ],
    hdrs = [
        "tsubstructure.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

pybind_extension(
    name = "lillymol_tools",
    srcs = [
        "lillymol_pybind_tools.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:dicer_api",
        "//Molecule_Tools:nvrtspsa",
        "//Molecule_Tools:unique_molecules_api",
        "//Molecule_Tools_Bdb:iwecfp_database_lookup_lib_shared",
        "//Molecule_Tools_Bdb:selimsteg",
        "//Molecule_Tools_Bdb:structure_database_shared",
        "@third_party//:berkeley_shared",
        "@python",
    ],
)

# Cannot figure out how to make this work. TODO:ianwatson
#py_test(
#    name = "test_lillymol",
#    size = "small",
#    srcs = ["lillymol_test.py"],
#    deps = [
#      ":lillymol_pybind",
#    ],
#)
