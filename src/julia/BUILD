cc_library (
    name = "lillymol_julia",
    srcs = [
        "lillymol.cc",
    ],
    deps = [
        ":resizable_array_holder",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:ring_data",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Tools:xlogp_lib",
        "//Molecule_Tools:alogp_lib",
        "@com_gitlab_libeigen_eigen//:eigen",
        "@libcxxwrap//:libcxxwrap",
        "@julia//:julia",
    ],
    tags = [
        "julia",
    ],
)

cc_library (
    name = "lillymol_substructure",
    srcs = [
        "lillymol_substructure.cc",
    ],
    deps = [
        ":resizable_array_holder",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@libcxxwrap//:libcxxwrap",
        "@julia//:julia",
    ],
    tags = [
        "julia",
    ],
)


cc_library(
    name = "bag_julia",
    srcs = [
        "bag.cc",
    ],
    deps = [
        "@libcxxwrap//:libcxxwrap",
        "@julia//:julia",
    ],
    tags = [
        "julia",
    ],
)

cc_shared_library(
    name = "lillymol_shared",
    deps = [
        ":lillymol_julia",
    ],
    dynamic_deps = [
    ],
    shared_lib_name = "lillymol_julia.so",
    tags = [
        "julia",
    ],
)

cc_shared_library(
    name = "lillymol_substructure_shared",
    deps = [
        ":lillymol_julia",
    ],
    dynamic_deps = [
    ],
    shared_lib_name = "lillymol_substructure_julia.so",
    tags = [
        "julia",
    ],
)


cc_shared_library(
    name = "box_shared",
    deps = [
        ":bag_julia",
    ],
    dynamic_deps = [
    ],
    shared_lib_name = "bag_julia.so",
    tags = [
        "julia",
    ],
)

cc_library (
    name = "resizable_array_holder",
    hdrs = [
        "resizable_array_holder.h"
    ],
    deps = [
      "//Foundational/iwaray",
    ],
    tags = [
        "julia",
    ],
)
