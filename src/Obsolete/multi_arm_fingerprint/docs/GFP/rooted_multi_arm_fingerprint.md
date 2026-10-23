# Rooted multi-arm fingerprint: design proposal

This document records the proposed design for a molecular fingerprint built from
three or four paths extending from a common central atom. An initial C++ generator and CLI are now implemented;
see [multi_arm_fingerprint](multi_arm_fingerprint.md) for current options and
defaults. The final fingerprint name and recommended defaults remain provisional.

“Rooted multi-arm fingerprint” is the working descriptive name. “Y-shaped
fingerprint” remains a useful shorthand, but does not describe the four-arm case.
“Star fingerprint” is another possible name, although the arms contain paths,
not just immediate neighbours. These shapes describe connectivity, not geometry.

## Purpose

The fingerprint captures combinations of paths around branching atoms. Unlike a
complete circular environment, a feature can describe selected extensions without
including every neighbouring branch along those paths. This may allow useful
partial matches between molecules whose complete local environments differ.

There is no assumption that it will outperform an extended-connectivity (EC)
fingerprint. The first goal is a clear, deterministic definition that can be
implemented efficiently and evaluated both alone and alongside existing LillyMol
fingerprints. This revisits the motivation behind the older `iwpathd` tool,
without adopting its enumeration or canonicalisation machinery.

## Centres and arms

* A centre has three or four heavy-atom neighbours.
* A degree-three centre generates three-arm features.
* A degree-four centre generates four-arm features, rather than enumerating its
  four possible three-arm subsets. This avoids an automatic combinatorial boost
  in contributions from degree-four centres.
* Other centre degrees are excluded initially.
* Each physical arm is identified by one distinct first neighbour of the centre.
* Every emitted feature includes one nonempty path from every physical arm.
  Two-atom and linear-path degeneracies are therefore excluded.
* Paths continue through branching atoms by enumerating the possible extensions.
  They do not describe the entire neighbourhood of those atoms.

Degree here refers to heavy-atom connectivity, not bond-order sum. The generator
rejects molecules containing explicit hydrogen atoms, checking once per call
before enumeration rather than checking each atom during traversal.

## Selecting roots with queries

By default, process every eligible degree-three and degree-four atom. Support
LillyMol's usual `-s` inline SMARTS and `-q` query-file options to restrict the
starting atoms.

When queries are supplied, collect the first matched atom (query atom zero) from
every embedding of every query. Take the union of these atoms, then apply the
centre-degree restrictions. Process each selected eligible root once, even if
multiple embeddings or queries select it. Query selection restricts roots only;
it does not restrict which atoms their arms may traverse. A selected atom with an
ineligible degree does not generate features.

A molecule matching no specified queries is an error by default. With `-z i`,
it remains in the output with an empty fingerprint; verbose statistics count
these query nonmatches. A query match selecting no eligible roots is allowed,
produces an empty fingerprint, and is not counted as a query nonmatch.

## Independent paths and rings

Each path maintains its own visited atoms, including the centre. It does not
immediately backtrack along the bond it just traversed. Otherwise, every eligible
extension is explored.

If an extension reaches an atom already on that path, that variant terminates
with a virtual ring-closure token. The token incorporates the closure bond type
and the size of the cycle closed: the distance along the current path back to
the encountered atom, plus one. It is not an explicit atom. Other extensions from
the same endpoint remain eligible, including other ring-closing extensions.

Visited atoms are independent between arms. Two arms can traverse the same atoms,
or traverse a ring in opposite directions. Consequently, a feature describes a
collection of rooted paths; it need not correspond to a tree-shaped molecular
subgraph. There is no cross-arm ownership or traversal-order dependence.

Every ordinary path prefix is eligible. A leaf or ring closure prevents further
extension of that variant, but its shorter path remains usable while other arms
extend. There is no padding or “stuck” token. Enumeration ends when there are no
more valid combinations within the length limits.

## Length definition

Lengths count bond steps outward from the centre. A ring-closing step counts as
one step, even though its terminal token is virtual. The total feature length is
the sum of its arm lengths; overlaps between arms are counted independently.

Every arm has length at least one. The smallest features therefore have total
length three or four, depending on the centre degree. A maximum total length B
bounds an individual arm at B - (k - 1), where k is the number of arms. An optional
maximum individual arm length provides a tighter limit.

## Two-phase enumeration at each centre

### Phase one: enumerate arms

Enumerate all eligible paths for each first neighbour, up to the effective arm
length limit. Compute hashes incrementally and store them in length buckets:

```
arm[first_neighbour][length] -> sequence of uint32_t hashes
```

The graph traversal, visited-atom tracking, and backtracking occur only in this
phase. Storage is reused or released after processing the centre, rather than
retaining all paths for the whole molecule.

Initially preserve every path occurrence, including equal hashes. Equal hashes
usually represent equivalent descriptions, but can also be hash collisions.
Later, identical (length, hash) entries could be compressed into multiplicities
if the counting semantics are preserved.

### Phase two: form fingerprint features

Enumerate arm-length tuples satisfying the minimum and maximum total length and
individual arm limits. For each tuple, enumerate the Cartesian product of its
hash buckets, selecting one entry from each physical arm. Canonicalise each
selection and increment its sparse fingerprint bit.

Physical arms are selected once each. Do not enumerate permutations of the arms:
canonical sorting handles their order. Different path selections that produce
the same bit still contribute separate occurrences.

This separates graph traversal from feature formation and avoids rediscovering
an arm path for each combination with other arms. It does not eliminate the
potentially large Cartesian product; length limits remain essential.

## Atom types, bond types, and hashing

Use `Atom_Typing_Specification` to assign atomic invariants once per molecule.
Initially favour local atom types so that the fingerprint does not implicitly
require agreement of extended neighbourhoods outside the selected paths.

Use fixed, distinct numeric bond codes. Test aromaticity before testing single,
double, or triple bond type, so aromatic encoding does not depend on an
underlying single/double representation. Define an explicit policy for any
remaining bond types before implementation is finalised.

An arm hash describes the alternating sequence of bond codes and atom invariants
outward from the centre. It also incorporates its length. A ring-closing terminal
step incorporates its bond code, a distinct closure marker, and ring size instead
of an ordinary terminal atom. Token encoding must distinguish atom, bond, and
closure fields before hashing; accidental 32-bit collisions remain possible.

For each selected combination:

1. Sort its three or four arm hashes as unsigned values.
2. Combine a feature seed, central atom invariant, arm count, and the ordered arm
   hashes using deterministic integer mixing.
3. Use the resulting `uint32_t` as a bit identifier in
   `Sparse_Fingerprint_Creator`, incrementing its occurrence count.

The central atom invariant is included once in the final combination. Retain
repeated arm hashes: `[a, a, b]` and `[a, b, b]` are different descriptions.
Equal hashes require no tie-breaker or special coefficient rule. Prefer sequential
mixing to a simple weighted sum, which can introduce avoidable algebraic
collisions.

All arithmetic must use explicitly defined unsigned widths and wraparound.
Hashing must not depend on atom numbers, neighbour iteration order, pointers,
platform byte order, or `std::hash`. Sorting makes the arm combination canonical;
the overall representation is still a hashed description, not a collision-free
canonical molecular graph. Select and document the exact mixing constants and
algorithm before generating reference fingerprints.

## Proposed command-line controls

Exact option letters and defaults should follow inspection of neighbouring
LillyMol generators. Proposed behavioural controls are:

| Control | Meaning |
| --- | --- |
| Atom typing | Existing `Atom_Typing_Specification` syntax. |
| Root queries (`-s`, `-q`) | Inline SMARTS or query files; the first atom of each embedding selects a root. Without queries, consider all eligible centres. |
| Maximum total length | Maximum sum of the arm lengths. |
| Maximum individual arm length | Optional additional limit. |
| Minimum total length | Exclude small features; every arm must still be nonempty. |
| Centre degree | Three, four, or both; initially propose both. |
| Ring closures | Include terminal closure tokens, or omit ring-closing extensions. Ordinary prefixes remain eligible in either mode. |
| Counts or presence | Accumulate occurrences, or reduce each final bit count to one. Initially propose counts. |

Validate incompatible bounds. A maximum total length below four can legitimately
exclude degree-four features while retaining degree-three features. Reject a
configuration that cannot generate either enabled feature class.

Provide the usual generator options for input type, aromaticity, elements,
standardisation, largest-fragment selection, fingerprint tag, verbosity, and TDT
filter mode. Normal output should include `$SMI<`, `PCN<`, and a sparse fingerprint
record; filter mode inserts fingerprints into existing TDT records. The interface
should support eventual deployment through `gfp_make.pl`. Tag and executable names
are not yet selected.

Verbose diagnostics should report eligible centres, paths enumerated, feature
combinations generated, and distinct bits. Do not silently truncate enumeration
at a feature-count threshold: that could introduce atom-order dependence. Any
future resource limit must have explicit failure semantics.

## Implementation structure and checks

Prefer reusable fingerprint-generation logic with a thin CLI wrapper. Keep
per-molecule atom typing and per-centre enumeration state separate from output
formatting. Do not revive unrelated complexity from `iwpathd` merely for reuse.

Focused correctness tests should cover:

* Atom-renumbering and neighbour-order invariance, including symmetric centres.
* Degree-three and degree-four features and the exclusion of incomplete arms.
* Root selection with `-s` and `-q`, multiple embeddings and queries, duplicate
  root selection, ineligible selected roots, and empty selections. Verify that
  queries restrict roots without restricting arm traversal.
* Exact total and individual length boundaries, leaf termination, and branching.
* Closure at the centre and at another earlier path atom, with correct ring size.
* Independent overlapping arms and multiple closure alternatives.
* Aromatic bond encoding across supported equivalent input representations.
* Repeated arm hashes, occurrence counts, and presence-only output.
* Reproducible fixed hash examples across Linux and macOS.
* Normal TDT output, filter mode, and an eventual `gfp_make.pl` integration case.

The enumeration and hash definition should be fixed before performance tuning.
Measure arm storage and combination counts separately to identify the actual
bottleneck. Multiplicity compression is a possible optimisation, not an initial
change in feature semantics.

## Evaluation

The primary proposed evaluation is LillyMol fingerprint-kernel SVM modelling.
Also assess conventional tabular models. Compare EC alone, multi-arm alone,
EC plus multi-arm, and EC plus linear paths, with comparable atom typing and
validation splits. Use scaffold-separated evaluation where practical.

Compare counted and binary features, degree-three versus degree-four centres,
length limits, and ring-token behaviour. Examine fingerprint density and runtime,
not only predictive performance. When combining fingerprints, control their
relative influence so a change in feature count does not obscure their value.
A separately weighted kernel combination can help with this comparison.

Questions left for measurements include useful default lengths, whether four-arm
features add value, the effect of independently overlapping arms, and whether the
new fingerprint improves upon or complements EC. No predictive advantage is
assumed in advance.

## Initial implementation decisions

The initial tool is `multi_arm_fingerprint`, with sparse tag `NCMARM`, total
lengths three through six, compressed atomic-number typing `UST:Y`, both centre
degrees, ring tokens enabled, and occurrence counts. See the linked tool document
for the concrete option letters. Python bindings and `gfp_make.pl` integration
remain deferred.

Hash version one uses MurmurHash3-style 32-bit word mixing and avalanche, without
byte reinterpretation. Arm seed is `0x41524d31`, atom marker `0x41544f4d`, closure
marker `0x52494e47`, and feature seed `0x4d415231`. Bond codes are aromatic=1,
single=2, double=3, triple=4, other=5. The source defines the exact arithmetic;
fixed explicit-atom-type regression examples retain its numerical results.
Arm hashes are finalised after mixing length; feature hashes are finalised after
mixing the centre type, arm count, and sorted arm hashes. Other bond types share
one category. Molecules containing explicit hydrogens are rejected before enumeration.

## Synchronised expansion experiment

An optional `-S` mode now restricts phase-two selections to a common radius,
retaining each shorter terminal path independently. At least one selected path
must reach the current radius, so selections whose arms have all terminated are
not counted again. Terminal status is recorded during phase one; a ring closure
is terminal, as is a path with no permitted continuation. Reaching a configured
length limit does not create a new chemical terminal token.

The total-length bounds still apply to actual lengths. `-L` optionally caps the
radius; for example, `-S -L 3 -R 12` permits all three rounds for both centre
degrees. Alternative paths of equal length can still generate multiple features.
The original enumeration remains the default for comparison. Hashing is unchanged.
This experiment addresses the concern that numerous large features dominate
similarity despite the continued presence of smaller, more transferable features.
