load("@nanobind_bazel//:build_defs.bzl", "nanobind_extension")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_python//python:defs.bzl", "py_test")

package(default_visibility = ["//visibility:public"])

genrule(
    name = "lillymol_soabi",
    outs = ["lillymol.so.soabi"],
    cmd = "$(PYTHON3) -c \"import sysconfig; print(sysconfig.get_config_var('SOABI'))\" > $@",
    toolchains = ["@rules_python//python:current_py_toolchain"],
)

nanobind_extension(
    name = "lillymol",
    srcs = [
        "lillymol_nb.cc",
        "lillymol_nb_atom_bond.cc",
        "lillymol_nb_chirality.cc",
        "lillymol_nb_common.cc",
        "lillymol_nb_descriptors.cc",
        "lillymol_nb_fingerprint.cc",
        "lillymol_nb_internal.h",
        "lillymol_nb_io.cc",
        "lillymol_nb_molecule.cc",
        "lillymol_nb_reaction.cc",
        "lillymol_nb_set_of_atoms.cc",
        "lillymol_nb_standardise.cc",
        "lillymol_nb_substructure.cc",
        "lillymol_nb_tools.cc",
        "lillymol_nb_tsubstructure.cc",
    ],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Molecule_Lib:atom_pair_fingerprint",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:linear_fingerprint",
        "//Molecule_Lib:moleculeio",
        ":tsubstructure",
        "//Molecule_Tools:alogp_lib",
        "//Molecule_Tools:dicer_api",
        "//Molecule_Tools:ec_fingerprint_lib",
        "//Molecule_Tools:iwdescr_lib_shared",
        "//Molecule_Tools:jwcats_lib",
        "//Molecule_Tools:medchemwizard_lib",
        "//Molecule_Tools:mformula",
        "//Molecule_Tools:nvrtspsa",
        "//Molecule_Tools:qed_lib",
        "//Molecule_Tools:ring_replacement_lib",
        "//Molecule_Tools:unique_molecules_api",
        "//Molecule_Tools:xlogp_lib",
        "//Utilities/GFP_Tools:gfp_context",
        "//Utilities/GFP_Tools:truncated_distance_matrix",
    ],
)

cc_library(
    name = "tsubstructure",
    srcs = ["tsubstructure.cc"],
    hdrs = ["tsubstructure.h"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

py_test(
    name = "lillymol_doc_test",
    srcs = ["lillymol_nb_doc_test.py"],
    main = "lillymol_nb_doc_test.py",
    data = [":lillymol"],
)

py_test(
    name = "lillymol_test",
    srcs = ["lillymol_nb_test.py"],
    main = "lillymol_nb_test.py",
    data = [
        ":lillymol",
        "testdata/ring_replacement_connectivity.textproto",
        "testdata/rand10.standard.gfp",
    ],
    deps = [
        "//Molecule_Lib:substructure_py_proto",
        "//Utilities/GFP_Tools:nearneighbours_py_proto",
    ],
)

nanobind_extension(
    name = "lillymol_gfp_server",
    srcs = ["lillymol_nb_gfp_server.cc"],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = ["//Utilities/GFP_Tools:gfp_server_lib"],
)

# Optional because //Depict:coords_2d depends on the separately installed
# coordgen library. build_linux.sh selects this target when
# BUILD_SCHRODINGER_2D is set.
nanobind_extension(
    name = "lillymol_depict",
    srcs = ["lillymol_nb_depict.cc"],
    dynamic_deps = ["//Molecule_Lib:molecule_lib_shared"],
    deps = [
        "//Depict:coords_2d",
        "//Depict:wedge_2d",
        "//Molecule_Lib:iwmolecule",
    ],
    tags = ["manual"],
)

py_test(
    name = "lillymol_nb_depict_test",
    srcs = ["lillymol_nb_depict_test.py"],
    data = [
        ":lillymol",
        ":lillymol_depict",
    ],
    tags = ["manual"],
)

py_test(
    name = "gfp_http_server_test",
    srcs = ["gfp_http_server_test.py"],
    data = ["testdata/rand10.gfp"],
    deps = [
        ":lillymol_gfp_server",
        "//Utilities/GFP_Tools:gfp_http_server",
        "//Utilities/GFP_Tools:nn_request_py_proto",
    ],
)


nanobind_extension(
    name = "lillymol_bdb",
    srcs = ["lillymol_nb_bdb.cc", "lillymol_nb_ring_database.cc"],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Tools_Bdb:iwecfp_database_lookup_lib_shared",
        "//Molecule_Tools_Bdb:ring_system_lookup_shared",
        "//Molecule_Tools_Bdb:selimsteg",
        "//Molecule_Tools_Bdb:structure_database_shared",
        "//Molecule_Tools_Bdb:substituent_identification_lib_shared",
        "@third_party//:berkeley_shared",
    ],
)

py_test(
    name = "lillymol_bdb_test",
    srcs = ["lillymol_nb_bdb_test.py"],
    main = "lillymol_nb_bdb_test.py",
    data = [
        ":lillymol",
        ":lillymol_bdb",
        "//Molecule_Tools_Bdb:smi2rings_bdb",
        "//BerkeleyDB:iwbdb_load",
        "//Molecule_Tools_Bdb:structure_database_load",
        "//Molecule_Tools_Bdb:substituent_identification",
    ],
)
