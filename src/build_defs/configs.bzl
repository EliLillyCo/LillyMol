config_setting(
    name = "gcc_omp",
    values = {
        "copts": "-fopenmp",
        "linkopts": "-lgomp",
    },
)

