load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library", "cc_test")

# 2D depiction. Currently coordinate generation only; a renderer is intended
# to follow, which is why this is its own package rather than part of
# Molecule_Lib. Keeping it separate also means coordgen is linked only into the
# tools that actually depict, instead of into every binary that depends on
# //Molecule_Lib:iwmolecule.
#
# Everything here needs @third_party//:coordgen, which only exists once
# src/build_linux.sh has been run with BUILD_SCHRODINGER_2D set. That is why the
# build script only descends into this package under the same condition.

cc_library(
    name = "coords_2d",
    srcs = [
        "coords_2d.cc",
    ],
    hdrs = [
        "coords_2d.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "@third_party//:coordgen",
    ],
    visibility = ["//visibility:public"],
)

# Wedge/hash bond assignment. Needs only coordinates that are already present,
# so unlike coords_2d it does not depend on coordgen.
cc_library(
    name = "wedge_2d",
    srcs = [
        "wedge_2d.cc",
    ],
    hdrs = [
        "wedge_2d.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
    visibility = ["//visibility:public"],
)

cc_test(
    name = "wedge_2d_test",
    timeout = "short",
    srcs = [
        "wedge_2d_test.cc",
    ],
    deps = [
        ":wedge_2d",
        # Only for the layouts the tests wedge - wedge_2d itself does not need it.
        ":coords_2d",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "coords_2d_test",
    timeout = "short",
    srcs = [
        "coords_2d_test.cc",
    ],
    deps = [
        ":coords_2d",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@googletest//:gtest",
        "@googletest//:gtest_main",
    ],
)

cc_binary(
    name = "make_2d_coordinates",
    srcs = [
        "make_2d_coordinates.cc",
    ],
    deps = [
        ":coords_2d",
        ":wedge_2d",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)
