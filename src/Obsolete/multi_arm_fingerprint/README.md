# Multi-arm fingerprint experiment

Archived on 2026-10-22. The observed similarity results did not justify further
work at this time. The relatively small number of three- and four-connected root
atoms in typical molecules is one possible explanation; it has not been
established as the cause.

This archive preserves the generator, CLI, unit tests, CLI regression fixtures,
user documentation, and design proposal, including the synchronised-expansion
variant. Source and fixture contents were moved without modification.

Layout and former locations:

- `src/`: formerly `src/Molecule_Tools/multi_arm_fingerprint*`.
- `docs/GFP/`: formerly the corresponding files under repository `docs/GFP/`.
- `test/multi_arm_fingerprint/`: formerly under repository `test/`.
- `BUILD.fragment`: the three former targets in `src/Molecule_Tools/BUILD`.

There is intentionally no active Bazel package or regression registration here.
To revisit the experiment, restore the files to their former locations and add
the saved build fragment to `src/Molecule_Tools/BUILD`. The CLI still references
that directory's `fingerprint_writer.h`; the saved build definitions retain its
original dependency. The archived documentation describes the original paths
and commands.

This experiment was not integrated into `gfp_make.pl` or the Python bindings.
