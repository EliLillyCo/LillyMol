cc_library(
    name = "lillymol_julia",
    srcs = [
        "lillymol.cc",
    ],
    tags = [
        "julia",
    ],
    deps = [
        ":resizable_array_holder",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:ring_data",
        "//Molecule_Tools:alogp_lib",
        "//Molecule_Tools:xlogp_lib",
        "@com_gitlab_libeigen_eigen//:eigen",
        "@julia",
        "@libcxxwrap",
    ],
)

cc_library(
    name = "lillymol_substructure",
    srcs = [
        "lillymol_substructure.cc",
    ],
    tags = [
        "julia",
    ],
    deps = [
        ":resizable_array_holder",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@julia",
        "@libcxxwrap",
    ],
)

cc_library(
    name = "bag_julia",
    srcs = [
        "bag.cc",
    ],
    tags = [
        "julia",
    ],
    deps = [
        "@julia",
        "@libcxxwrap",
    ],
)

cc_shared_library(
    name = "lillymol_shared",
    dynamic_deps = [
    ],
    shared_lib_name = "lillymol_julia.so",
    tags = [
        "julia",
    ],
    deps = [
        ":lillymol_julia",
    ],
)

cc_shared_library(
    name = "lillymol_substructure_shared",
    dynamic_deps = [
    ],
    shared_lib_name = "lillymol_substructure_julia.so",
    tags = [
        "julia",
    ],
    deps = [
        ":lillymol_julia",
    ],
)

cc_shared_library(
    name = "box_shared",
    dynamic_deps = [
    ],
    shared_lib_name = "bag_julia.so",
    tags = [
        "julia",
    ],
    deps = [
        ":bag_julia",
    ],
)

cc_library(
    name = "resizable_array_holder",
    hdrs = [
        "resizable_array_holder.h",
    ],
    tags = [
        "julia",
    ],
    deps = [
        "//Foundational/iwaray",
    ],
)
