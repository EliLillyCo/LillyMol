# multi_arm_fingerprint

`multi_arm_fingerprint` generates sparse counted fingerprints from three or four
paths extending from a common atom. Every feature includes all three or four
arms. Degree-four centres produce four-arm features, rather than four subsets
of three arms. The generator is experimental; predictive value and useful length
defaults remain to be evaluated.

```
multi_arm_fingerprint -R 6 -P UST:Y molecules.smi > molecules.gfp
multi_arm_fingerprint -R 8 -L 3 -s '[ND3]' molecules.smi > nitrogen_roots.gfp
multi_arm_fingerprint -q PROTO:roots.textproto molecules.smi > selected.gfp
```

Default output contains `$SMI<`, `PCN<`, and `NCMARM<` records followed by `|`.
Molecules with no eligible roots remain in the output with `NCMARM<>`.
`-J` uses the common LillyMol fingerprint writer, including its fixed-width and
array output controls (`-J help`). `-f` reads TDT records and inserts fingerprints
after their `$SMI<` lines. `-` reads standard input; outside filter mode its default
input type is SMILES.

| Option | Meaning |
| --- | --- |
| `-S` | Synchronised expansion: combine paths at a common radius, retaining shorter terminal paths. Default is all length combinations. |
| `-R n` | Maximum total arm length, default 6. |
| `-r n` | Minimum total arm length, default 3. |
| `-L n` | Maximum individual arm length; default bounded by total length. |
| `-d 3` or `-d 4` | Restrict centre degree; default both. |
| `-c` | Omit ring-closing extensions; ordinary prefixes remain eligible. |
| `-b` | Set each final bit count to one. |
| `-P specification` | LillyMol atom typing; default `UST:Y` (compressed atomic number). |
| `-s SMARTS` | Select roots using the first atom of every query embedding. Repeatable. |
| `-q query` | Select roots using query files and standard LillyMol query-file syntax. Repeatable. |
| `-z i` | Allow molecules matching no queries and write empty fingerprints; otherwise a query nonmatch is an error. |
| `-J tag` | Fingerprint output specification, default `NCMARM`. |
| `-f` | TDT filter mode; requires fingerprint rather than descriptor/SVML output. |
| `-l` | Reduce to largest fragment before root selection and generation. |
| `-g specification` | Chemical standardisation before root selection and generation. |
| `-i type` | Input type; otherwise inferred from filename. |
| `-A specification` | Standard aromaticity options. |
| `-E specification` | Standard element options. |
| `-v` | Report enumeration totals and minimum, maximum, and average distinct bits per molecule (including zero-bit molecules). |

Lengths count bonds. A ring-closing step counts as one, including its terminal
virtual token. Every arm must be nonempty, so degree-three features require at
least three total steps and degree-four features require at least four. Degree
means heavy-atom connectivity. Molecules containing explicit hydrogen atoms are
rejected; remove explicit hydrogens before fingerprinting.
Queries select a union of roots, each processed once. With queries specified,
a molecule matching none is an error unless `-z i` is given. With `-z i`, retain
the record with an empty fingerprint; verbose output reports the number of
nonmatching molecules. A match selecting only ineligible roots is allowed and
is not counted as a query nonmatch. Degree restrictions still
apply. Queries do not restrict the atoms traversed by the arms.

Each arm independently tracks visited atoms. Immediate reversal of its last bond
is excluded. An extension to an earlier atom terminates that variant with a token
encoding the bond type and the size of the cycle just closed. Different arms may
traverse the same atoms. Aromatic bonds are encoded before single/double bonds;
stereochemistry is not explicitly encoded.

The implementation first enumerates paths for each arm into length buckets, then
combines one path from each arm within the total-length budget. Sorted arm hashes
are combined with the centre atom type and number of arms. Distinct occurrences
of the same description increment the same bit, unless `-b` is specified. Sparse
fingerprint serialisation uses the existing LillyMol count representation; its
byte-sized counts saturate at 255. Large length limits can generate many paths
and combinations; no combinations are silently discarded to control runtime.

## Synchronised expansion

Use `-S` to reduce the combinations of different arm lengths:

```
multi_arm_fingerprint -S -L 3 -R 12 molecules.smi > synchronised.gfp
```

At radius one, select a one-bond path from every arm. At each later radius,
select paths of that length, or shorter paths that have terminated at a leaf or
ring closure. Each terminal variant remains eligible independently, even when
other variants of the same arm continue. A combination is emitted only at the
radius of its longest path, preventing repetition once all its arms have stopped.

`-R` still bounds the sum of actual path lengths, and `-r` still bounds the
minimum sum. `-L` caps the radius. To allow every four-arm combination through
radius three, use `-L 3 -R 12`; with `-R 6`, four arms of length two would exceed
the total budget. If `-L` is omitted, the total budget bounds the radius.

Alternative paths within an arm still produce Cartesian combinations at each
radius. This mode removes arbitrary combinations between lengths, rather than
guaranteeing one bit per root per radius. For example, a centre with three
unbranched two-bond arms generates two combinations instead of the original
eight when the total budget permits both radii. `-b` remains available to reduce
final bit counts to presence only. Without `-S`, enumeration and bit identities
remain unchanged. Identical path selections share the same bit identities across
both modes, enabling direct comparisons.

The C++ control is `Options::_synchronised`, default `false`. The generator still
uses two phases, recording terminal status alongside each arm hash in phase one.
Whether this improves molecular similarity or predictive performance remains to
be evaluated.

## C++ interface

`Molecule_Tools/multi_arm_fingerprint.h` provides
`multi_arm_fingerprint::Generator`, independent of command-line processing and
output. It supports configurable atom typing and an overload taking caller-supplied
`uint32_t` atom invariants. An optional integer mask restricts root atoms; optional
`Statistics` reports the work performed for that molecule. Each successful call
replaces the destination `Sparse_Fingerprint_Creator`, including for empty results.

```
multi_arm_fingerprint::Generator generator;
multi_arm_fingerprint::Options options;
options._max_total_length = 6;
if (!generator.SetOptions(options)) {
  // Invalid configuration.
}
Sparse_Fingerprint_Creator fingerprint;
if (!generator.Fingerprint(molecule, fingerprint)) {
  // Generation failed; fingerprint is empty.
}
```

The public Bazel library target is
`//Molecule_Tools:multi_arm_fingerprint_lib`. Python bindings are deferred until
the CLI and fingerprint have been evaluated; the same generator is intended to
serve those bindings. The tool is not yet integrated into `gfp_make.pl`.

## Building and testing

From `src`:

```
bazel build -c opt --config=native //Molecule_Tools:multi_arm_fingerprint
bazel test -c opt --config=native //Molecule_Tools:multi_arm_fingerprint_test
ruby ../test/run_all_test.rb -bindir "$PWD/bazel-bin/Molecule_Tools" multi_arm_fingerprint
```

See [the design document](rooted_multi_arm_fingerprint.md) for the rationale and
proposed modelling comparisons.
