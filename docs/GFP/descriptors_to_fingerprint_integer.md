# Descriptors To Fingerprint Integer

`descriptors_to_fingerprint_integer` converts a tabular file of integer descriptors into a LillyMol GFP/TDT fingerprint stream.

It is intended for cases where descriptor values are already available
as integer vector components and there is a need to incorporate this
information into a similarity calculation.
The generated GFP/TDT output can be used directly with tools such as
gfp_lnearneighbours, gfp_pairwise_distances, gfp_spread, and other
LillyMol fingerprint processing utilities.

The sparse fingerprint form preserves
integer counts and is suitable for count-based similarities such as
Tanimoto similarity on integer vectors.

## Input format

The descriptor file is a whitespace-separated table by default.

The first row is a header. The first column is the identifier column, followed by one column per descriptor.

Example:

```text
Name F1 F2 F3 F4
active1 1 0 3 4
active2 2 0 2 5
inactive1 0 7 1 0
```

All descriptor values must be non-negative integers in the range 0 to 255.

By default, fields are separated by spaces. Use `-i tab` for tab-separated input.

## Basic usage

```bash
descriptors_to_fingerprint_integer -s C descriptors.dat > descriptors.gfp
```

This writes a GFP/TDT file with a dummy SMILES value for each record.

Example output:

```text
$SMI<C>
PCN<active1>
NCDSCI<...>
|
```

Important:  Sparse fingerprint bit assignments are generated
dynamically as non-zero descriptor values are encountered.  If
multiple descriptor files will later be compared, they must be
processed in a single invocation so that descriptor columns are
assigned consistent bit numbers.

Given identical input files in the same order, bit assignments are deterministic and reproducible.

## Supplying SMILES

A SMILES file can be supplied with `-S`.

```bash
descriptors_to_fingerprint_integer -S molecules.smi descriptors.dat > descriptors.gfp
```

The SMILES file should contain SMILES and identifier fields:

```text
CC methane
CCC propane
```

The identifiers in the SMILES file must match the identifiers in the descriptor table.

## Sparse Counted Fingerprints

By default the tool writes sparse count fingerprints using the tag `NCDSCI<`.

```bash
descriptors_to_fingerprint_integer descriptors.dat > descriptors.gfp
```

Sparse fingerprints preserve the integer descriptor values. Zero-valued features are not written.

This is the preferred mode when the descriptor values are counts or
integer-valued features and will be compared with count-based
fingerprint similarity methods.

The tag can be changed with `-J`. 

```bash
descriptors_to_fingerprint_integer -J NCMYFP descriptors.dat > descriptors.gfp
```

Tags beginning with `NC` generate sparse, non-colliding fingerprints.

## Bit replication

For sparse fingerprints, `-r` replicates each descriptor column into multiple adjacent sparse bits.
This is useful with SVMFP models where there is a need to give increased
weight to one or more features. By default, each non zero feature sets
one bit, with the count associated with the bit being the value of the feature.
If the `-r` option is specified, that number of separate bits are generated.
So if the feature in column 23 maps to bit number 6 and -r 4 is in effect
a non zero value in column 23 might map to bits 6, 7, 8 and 9. The
influence of this feature in the similarity computation is multiplied by 4.

```bash
descriptors_to_fingerprint_integer -r 4 descriptors.dat > descriptors.gfp
```

If descriptor column `i` has value `v`, then `-r 4` writes four sparse
features derived from that column, each with count `v`.

This option only applies to sparse fingerprints.

## Fixed-width fingerprints

If the tag supplied with `-J` starts with `FP`, the tool writes a fixed-width fingerprint.

```bash
descriptors_to_fingerprint_integer -J FPDESC descriptors.dat > descriptors.gfp
```

In this mode, descriptor counts are not preserved.  Any non-zero
descriptor value sets the corresponding bit in the binary fingerprint.
Zero descriptor values leave the corresponding bit unset.

The number of bits is rounded up to a multiple of 8.

## TDT filter mode

The tool can also insert descriptor fingerprints into an existing GFP/TDT stream.

```bash
  generate fp | \
           descriptors_to_fingerprint_integer -D descriptors.dat -f - > augmented.gfp
```

In this mode:

* `existing.gfp` is read as a TDT stream.
* Each record is copied to output.
* When a `PCN<identifier>` record is encountered, the matching descriptor vector from `descriptors.dat` is inserted as a fingerprint record.

The descriptor file supplied with `-D` must contain identifiers that match the `PCN` identifiers in the input GFP/TDT file.

given an input GFP that might look like
```text
$SMI<c1ccccc1>
PCN<benzene>
FPXX<.....>
|
```
and with a command like
```bash
descriptors_to_fingerprint_integer -D descriptors.dat -f existing.gfp > augmented.gfp
```
Might generate output like
```text
$SMI<c1ccccc1>
PCN<benzene>
FPXX<.....>
NCDSC<....>
|
```

## Options

```text
-f             Work as a TDT filter, inserting fingerprints into an existing GFP/TDT stream.
-D <fname>     Descriptor file used in filter mode.
-S <fname>     SMILES file. Identifiers are matched against descriptor identifiers.
-s <smiles>    Use a constant dummy SMILES for each output record.
-J <tag>       Fingerprint tag. Tags starting with NC produce sparse count fingerprints.
               Tags starting with FP produce fixed-width fingerprints.
-r <n>         Number of sparse-bit replicates per descriptor column.
-i <char>      Input token separator. Use '-i tab' for tab-separated input.
-W <dir>       When processing multiple input files, generate separate output .gfp files with the same prefixes in <dir>.
-v             Verbose output.
```

## Example

Descriptor file:

```text
Name F1 F2 F3
active1 1 2 0
active2 1 3 0
inactive1 0 0 7
```

Command:

```bash
descriptors_to_fingerprint_integer -s '*' descriptors.dat > descriptors.gfp
```

The resulting GFP/TDT file can be used with LillyMol GFP similarity tools.

## Notes

The sparse fingerprint representation is especially useful for integer vectors where the magnitude of each feature matters. For example, two vectors can be compared using a count-based Tanimoto similarity where the numerator is the sum of element-wise minima and the denominator is the total count in both vectors minus that common count.

Fixed-width fingerprints should be used only when presence/absence information is desired.

There are some example files in [examples](descriptors_to_fingerprint_integer_example/README.md)
