# load("@rules_python_pytest//python_pytest:defs.bzl", "py_pytest_test")
# load("@rules_python//python:py_test.bzl", "py_test")
load("@pybind11_bazel//:build_defs.bzl", "pybind_extension")
load("@rules_python//python:defs.bzl", "py_library", "py_test")
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_proto_library", "cc_library", "cc_test")
load("@com_google_protobuf//bazel:cc_proto_library.bzl", "cc_proto_library")

package(default_visibility = ["//visibility:public"])

exports_files(["testdata/rand10.standard.gfp"])

pybind_extension(
    name = "lillymol",
    srcs = [
        "lillymol_pybind.cc",
    ],
    copts = [
        "-save-temps",
    ],
    # Molecule_Lib comes from the one shared copy, not statically per module.
    # See the comment on //Molecule_Lib:molecule_lib_shared.
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        ":molecule",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:alogp_lib",
        "//Molecule_Tools:xlogp_lib",
#       "@python",
#       "@third_party//:python",
    ],
)

pybind_extension(
    name = "lillymol_atom",
    srcs = [
        "lillymol_pybind_atom.cc",
    ],
    deps = [
        ":lillymol",
    ],
)

pybind_extension(
    name = "lillymol_bond",
    srcs = [
        "lillymol_pybind_bond.cc",
    ],
    deps = [
        ":lillymol",
    ],
)

pybind_extension(
    name = "lillymol_fingerprint",
    srcs = [
        "lillymol_pybind_fingerprint.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        ":molecule",
        "//Molecule_Lib:atom_pair_fingerprint",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:linear_fingerprint",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:ec_fingerprint_lib",
#       "@python",
    ],
)


pybind_extension(
    name = "lillymol_query",
    srcs = [
        "lillymol_pybind_query.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        ":molecule",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
#       "@python",
    ],
)

pybind_extension(
    name = "lillymol_standardise",
    srcs = [
        "lillymol_pybind_std.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
#       "@python",
    ],
)

pybind_extension(
    name = "lillymol_set_of_atoms",
    srcs = [
        "lillymol_pybind_set_of_atoms.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
#       "@python",
    ],
)

pybind_extension(
    name = "lillymol_ring",
    srcs = [
        "lillymol_pybind_ring.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
#       "@python",
    ],
)

pybind_extension(
    name = "lillymol_io",
    srcs = [
        "lillymol_pybind_io.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
#       "@python",
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
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
#       "@python",
    ],
)

pybind_extension(
    name = "lillymol_tsubstructure",
    srcs = [
        "lillymol_pybind_tsubstructure.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
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
    name = "lillymol_gfp_server",
    srcs = [
        "lillymol_pybind_gfp_server.cc",
    ],
    # Molecule_Lib arrives transitively through gfp_server_lib.
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Utilities/GFP_Tools:gfp_server_lib",
    ],
)

pybind_extension(
    name = "lillymol_tools",
    srcs = [
        "lillymol_pybind_tools.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools:alogp_lib",
        "//Molecule_Tools:dicer_api",
        "//Molecule_Tools:iwdescr_lib_shared",
        "//Molecule_Tools:jwcats_lib",
        "//Molecule_Tools:medchemwizard_lib",
        "//Molecule_Tools:nvrtspsa",
        "//Molecule_Tools:qed_lib",
        "//Molecule_Tools:ring_replacement_lib",
        "//Molecule_Tools:unique_molecules_api",
        "//Molecule_Tools:xlogp_lib",
        "//Molecule_Tools_Bdb:iwecfp_database_lookup_lib_shared",
        "//Molecule_Tools_Bdb:selimsteg",
        "//Molecule_Tools_Bdb:structure_database_shared",
        "//Utilities/GFP_Tools:gfp_context",
        "//Utilities/GFP_Tools:truncated_distance_matrix",
        "@third_party//:berkeley_shared",
#       "@python",
    ],
)

py_test(
    name = "lillymol_amw_test",
    srcs = ["lillymol_amw_test.py"],
    deps = [
        ":lillymol",
    ],
)

py_test(
    name = "lillymol_substructure_test",
    srcs = ["lillymol_substructure_test.py"],
    deps = [
        ":lillymol",
        ":lillymol_query",
        "//Molecule_Lib:substructure_py_proto",
    ],
)

py_test(
    name = "lillymol_tools_gfp_test",
    srcs = ["lillymol_tools_gfp_test.py"],
    data = [
        "testdata/rand10.complex.gfp",
        "testdata/rand10.standard.gfp",
        "//Utilities/GFP_Tools:gfp_nearneighbours_single_file_tbb",
    ],
    deps = [
        ":lillymol",
        ":lillymol_tools",
        "//Utilities/GFP_Tools:nearneighbours_py_proto",
    ],
)

py_test(
    name = "lillymol_tools_iwdescr_test",
    srcs = ["lillymol_tools_iwdescr_test.py"],
    deps = [
        ":lillymol",
        ":lillymol_tools",
    ],
)

py_test(
    name = "lillymol_tools_qed_test",
    srcs = ["lillymol_tools_qed_test.py"],
    deps = [
        ":lillymol",
        ":lillymol_tools",
    ],
)

py_test(
    name = "lillymol_tools_medchemwizard_test",
    srcs = ["lillymol_tools_medchemwizard_test.py"],
    deps = [
        ":lillymol",
        ":lillymol_tools",
    ],
)

py_test(
    name = "gfp_http_server_test",
    srcs = ["gfp_http_server_test.py"],
    data = [
        "testdata/rand10.gfp",
    ],
    deps = [
        ":lillymol_gfp_server",
        "//Utilities/GFP_Tools:gfp_http_server",
        "//Utilities/GFP_Tools:nn_request_py_proto",
    ],
)

py_test(
    name = "gfp_http_client_test",
    srcs = ["gfp_http_client_test.py"],
    deps = [
        "//Utilities/GFP_Tools:gfp_http_client",
        "//Utilities/GFP_Tools:nn_request_py_proto",
    ],
)

# Cannot figure out how to make this work. TODO:ianwatson
# Observe that Codex ran the test above by including
# --strategy=TestRunner=local --test_env=LILLYMOL_HOME=/home/ian/LillyMol_IWFORK
#py_test(
#    name = "test_lillymol",
#    size = "small",
#    srcs = ["lillymol_test.py"],
#    deps = [
#      ":lillymol_pybind",
#    ],
#)
