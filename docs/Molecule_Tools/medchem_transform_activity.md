# medchem_transform_activity

`medchem_transform_activity.rb` generates medchem-wizard products, looks them up
in a structure database, and summarizes activity changes for each reaction.
Measured activities are used for products found in the database; model predictions
can supply activities for products not found there.

## Generate, score externally, and restart

When scoring is expensive, run generation and database lookup separately:

```shell
medchem_transform_activity.rb generate -d structures.bdb input.smi
```

This creates the starting-molecule database and leaves `found.smi` and
`notfound.smi` in the current directory. It exits after lookup without scoring or
analysing the products. Use `--found` and `--notfound` to choose other paths;
`.smi` is appended if the supplied path does not already end with that suffix.
Keep each run's intermediate files in a separate directory or use distinct names.

Score `notfound.smi` using your preferred model or batch-processing system.
Scoring can run in parallel batches. Combine the results into `notfound.pred`
with exactly one header, preserving every product's full identifier:

```text
id prediction
mol1__MCTA__add_oxygen 3.0
mol2__MCTA__add_nitrogen 1.5
```

Prediction files require a header. Values are predicted absolute activities,
using the same scale as the measured activity file. The combined identifier
contains the starting ID and reaction name separated by `__MCTA__`; do not
truncate it to the starting ID. Reaction names cannot contain whitespace, and
neither component may contain the reserved separator.

Resume analysis once predictions are ready:

```shell
medchem_transform_activity.rb restart -A activity.txt \
  --found found.smi --predictions notfound.pred -o results.txt
```

`restart` is an alias for `analyse`. It reads existing found structures and
predictions without rebuilding the database, regenerating molecules, looking up
products, or invoking a model. It does not require the original SMILES file or
`-d`. Omitting `--predictions` produces measured-activity summaries only.

The original SMILES input has two columns, `smiles id`, without a header. The
activity file has a header followed by `id activity`; `.csv` files are parsed as
CSV and other files as whitespace-separated text. Positive reported changes mean
the transformed molecule has a higher numerical activity than its starting
molecule. Whether that is desirable depends on the activity scale.

## Run all stages together

```shell
medchem_transform_activity.rb run -A activity.txt -d structures.bdb \
  --model-command 'score_model {input} > {output}' input.smi
```

`run` is also the default mode. The command's `{input}` and `{output}` placeholders
are replaced with shell-escaped paths. Explicitly named intermediate files are
retained; unnamed temporary files are normally removed. Use `--keep` to retain
temporary files.

Use `--profile-json profiles.json --reactions /path/to/REACTIONS` with `run`,
`analyse`, or `restart` to additionally write reaction profiles containing
measured and predicted effects.

For an offline report with molecular examples, distributions, sorting, and count
filters, see [Medchem transformation activity report](medchem_transform_activity_html.md).

## Sorting candidates for prediction

After database lookup, `generate` and `run` sort the not-found molecules by
ascending atom count using `msort_parallel -k natoms`. The sorted structures
replace the contents of the same not-found file, so external scoring and
`--model-command` consume the sorted candidates. Combined molecule/reaction
identifiers are retained. All candidates are retained unless per-reaction sampling
is enabled as described below.

Sorting uses eight threads by default. Override the thread count or executable:

```shell
medchem_transform_activity.rb generate -d structures.bdb \
  --sort-threads 4 --msort-parallel /path/to/msort_parallel input.smi
```

For small files, the thread count is capped at the number of nonblank records.
Empty files are left untouched. Sorting writes to a temporary file in the same
directory and replaces the not-found file only after success and verification
that its byte size matches the original. A sorting failure or size mismatch
stops processing and leaves the original lookup output intact. `analyse` and
`restart` do not invoke the sorter.

## Sampling candidates by reaction

Use `--max-per-reaction N` to retain at most N products for each reaction from
the sorted not-found file:

```shell
medchem_transform_activity.rb generate -d structures.bdb \
  --max-per-reaction 100 --sort-threads 8 input.smi
```

Sampling runs after atom-count sorting in both `generate` and `run`, so the
smallest products are retained first. This favours simpler molecules for estimating
transformation effects. It does not limit reuse of a starting molecule.
No sampling is performed unless the option is supplied; N must be positive.

The script invokes `unique_rows -c 2 -r N -subset '.*__MCTA__(.*)'` to count only
the reaction portion of each combined identifier. Different starting identifiers
therefore share the same per-reaction allowance. Full identifiers, structures,
and sorted order remain intact. `--unique-rows EXE` selects a different executable;
otherwise `unique_rows` is resolved through PATH.

The sampled output replaces the same not-found file only after the command
succeeds. Unlike sorting, sampling intentionally reduces file size. If sampling
fails, the sorted file remains intact and scoring does not start. Empty files
skip sampling. Found molecules and their experimental summaries are unaffected;
`analyse` and `restart` do not sample existing results.
