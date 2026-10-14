# mispredicted

`mispredicted` combines prediction records from multiple files and identifies
molecules with large average prediction errors. These can be candidates for
investigating activity cliffs or measurement problems. The tool reports prediction
performance; it does not determine why a molecule is poorly predicted.

```shell
mispredicted -E measured.txt -n 2 -w 100 -d pred1.txt pred2.txt pred3.txt
```

This writes at most 100 identifiers with at least two valid prediction records,
ranked by mean absolute error. Identifiers below the minimum are excluded before
ranking and selection.

## Input files

Experimental and prediction files require a header. The first column is the
identifier; the second is the measurement or prediction by default:

```text
ID measured
mol1 2.0
mol2 4.0
```

```text
ID predicted
mol1 3.0
mol2 3.5
```

Use `-e <column>` and `-p <column>` to select other experimental and prediction
columns. Column numbers are one-based and must be at least 2. Continuous input
should use spaces as separators; classification input also translates tabs.
Identifiers are matched exactly and must be single tokens. Measurements and
predictions must use the same activity scale.

Repeated experimental identifiers are errors unless `-h` is specified; `-h`
retains the first value. Prediction files may contain different subsets of the
experimental identifiers. `-g` ignores predictions with no experimental match.

`-F` skips invalid prediction records and continues reading the file, retaining
valid records before and after the error. Invalid records do not contribute to
prediction counts. It does not relax experimental-data validation.

## Error metrics and selection

Continuous errors default to **mean absolute error**:

```text
mean(abs(measured - predicted))
```

`-w <count>` selects at most that many molecules with the largest mean errors.
`-W <threshold>` retains molecules whose mean error magnitude is at least the
threshold. Use one of these selection options at a time; if both are given,
`-w` takes precedence.

`-s` instead uses **mean signed error**, `mean(measured - predicted)`. Positive
values indicate under-prediction; negative values indicate over-prediction.
With `-s`, both `-w` and `-W` select by the magnitude of that signed mean.
Opposite-sign errors can cancel: errors of +2 and -2 have mean absolute error 2
but mean signed error 0.

`-b` requires `-w` and implies signed errors. It selects from both directions,
approximately balancing accumulated error magnitudes while both signs remain
available, and never exceeds the requested count. One-sided data remains valid.

`-n <count>` requires at least that many **valid prediction records** per
identifier, default 1. Duplicate records within a file count separately: this is
not a count of independent files or models. A filter with no qualifying molecules
succeeds with no data rows. Input containing no valid predictions is an error.

Continuous output sorts by decreasing mean error, then decreasing prediction
count, then ascending identifier. For signed errors this puts positive means
before negative means; selection still considers both directions. `-r` reverses
the final output order, including ties, without changing selected identifiers.

## Output

Default continuous output is a descriptive line per identifier. `-d` produces a
headered table:

```text
ID EXPT AVDIFF
mol1 2 1
mol2 4 0.5
```

`AVDIFF` is the selected mean error metric. `-S molecules.smi` prepends SMILES to
continuous output; the SMILES file contains `smiles id` records without a header.
Missing SMILES are reported and represented by `*`.

`-Y summary.txt` writes a continuous summary with these columns:

| Column | Meaning |
| --- | --- |
| `ID` | Identifier |
| `Expt` | Experimental value |
| `N` | Number of valid prediction records |
| `minpred`, `maxpred`, `avepred` | Prediction range and mean |
| `avediff` | Mean error using the selected metric |
| `maxdiff` | Largest error using the selected metric; in signed mode, the largest signed value |
| `stddiff` | Sample standard deviation of those errors |

A single prediction has `stddiff` 0 by convention. With the default metric this
is spread in absolute errors; with `-s` it is spread in signed errors.

`-T results` writes one single-line `MispredictedContinuous` textproto per
identifier to `results.textproto`. The suffix is appended if needed. Summary and
textproto outputs use the same eligible identifiers and order as primary output.
Output paths that alias input files or each other are rejected before opening.

## Selecting prediction files

Pass filenames explicitly, or use `-R` to match filenames in the current directory:

```shell
mispredicted -E measured.txt -n 5 -d -R '^prediction_[0-9]+\.txt$'
```

Matching is against the filename, not a recursive directory search. With `-R`,
positional prediction filenames are not used.

Continuous mode also accepts a list file:

```shell
mispredicted -E measured.txt -d F:prediction_files.txt
```

The list contains one filename per line, without a header. Empty lines and lines
starting with `#` are skipped. Paths are relative to the working directory.

## Classification

`-C` compares class-label tokens instead of numeric errors:

```shell
mispredicted -C -E measured_classes.txt -n 2 -d pred1.txt pred2.txt
```

Descriptor output has `ID EXPT N Correct Accuracy`, where accuracy is the fraction
of prediction records matching the experimental class. Default ordering is best
accuracy first, then largest prediction count, class label (length before lexical
order), and identifier. `-r`
reverses that order, putting the least accurate identifiers first.

`-c <cutoff>` enables classification from numeric scores: scores below the cutoff
become `-1`, and scores at or above it become `1`. `-X <label1,label2>` translates
predicted labels to experimental classes by prevalence: the first predicted label
maps to the more common experimental class. This requires exactly two experimental
classes; equal class prevalence currently has ambiguous mapping.

Classification currently does not implement `-w`, `-W`, `-Y`, or `-T`. Do not use
`-S` with classification descriptor output: it adds a SMILES header without the
corresponding row values. These limitations remain tracked in
[mispredicted_review.md](mispredicted_review.md).

## Support-vector distance options

`-k` reads and validates an additional distance column immediately after the
selected prediction column, in the range 0 to 1. Distances currently do not affect
ranking or appear in output. The advertised `-D` distance/error output option is
not implemented.

## Build and regression tests

From `src`:

```shell
bazel build -c opt --config=native //Utilities/General:mispredicted
ruby ../test/run_all_test.rb -bindir "$PWD/bazel-bin/Utilities/General" mispredicted
```

The native configuration is for local builds targeting the host CPU.
