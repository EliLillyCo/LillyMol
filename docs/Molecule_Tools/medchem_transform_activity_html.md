# Medchem transformation activity report

`medchem_transform_activity_html.sh` creates a self-contained HTML report for
exploring which chemical transformations increase or decrease a numerical response.
It is separate from `medchem_transform_activity.rb` and does not change its files.

```shell
medchem_transform_activity_html.sh \
  --starting-smiles input.smi --activity activity.txt \
  --found found.smi --notfound notfound.smi --predictions notfound.pred \
  --max-examples 5 -o transformations.html
```

Open `transformations.html` in a modern browser. All data, drawings, styles, and
controls are embedded: no server, external network requests, or browser-side
chemistry library is needed. The HTML can be shared as a single file.

## Evidence and structures

The generator reads the underlying records, rather than the summary table,
because distributions and molecular examples require individual observations.

- Starting SMILES: `smiles id` records, without a header.
- Activity: a header followed by `id activity`. Files ending in `.csv` use CSV;
  otherwise fields are separated by whitespace.
- Found products: `smiles id__MCTA__reaction matched_id`, without a header.
  Legacy `smiles id reaction matched_id` records are also accepted.
- Notfound products: `smiles id__MCTA__reaction`, without a header.
  Legacy `smiles id reaction` records are also accepted.
- Predictions: a header followed by `id__MCTA__reaction value`.
  Legacy `id reaction value` records are also accepted.

Experimental changes are matched-product measured activity minus starting measured
activity. Predicted changes are product predicted activity minus starting measured
activity. A positive change means a larger numerical response; whether that is
desirable depends on the assay. These two evidence sources remain visibly separate.

Either experimental or predicted evidence can be supplied alone. Predicted evidence
requires both `--notfound` and `--predictions`, because predictions contain no
structures. Predictions must match an existing starting ID and notfound product.
Notfound records without predictions are excluded and counted in the report.
Repeated prediction records all contribute; the number of distinct starting
molecules is shown separately from record count.

LillyMol performs SMILES parsing, atom counting, and structural identity checks.
RDKit is used only to draw the retained examples. Generation requires the staged
LillyMol Python bindings and RDKit in the same Python environment. The shell
wrapper invokes LillyMol's `run_python.sh` to check the Python ABI and locate the
staged libraries. Set `PYTHON` if a specific compatible interpreter is needed.

## Smallest examples and sorted input

The main pair is the smallest product exemplar across both evidence sources.
Additional pairs are also selected by product atom count, rather than activity
change. Starting structures are recovered by ID; they need not be sorted.

The default cap is five distinct pairs per transformation, including the main
pair. Set `--max-examples` to another positive number. Expanding a row renders
only its retained examples; it does not reveal an unlimited list. The cap does
not reduce the records used for statistics or distributions.

Input order need not be sorted by default. For large files, sorting products with
LillyMol first allows a faster path that avoids unnecessary structure processing:

```shell
msort -k natoms found.smi > found_sorted.smi
msort -k natoms notfound.smi > notfound_sorted.smi
medchem_transform_activity_html.sh \
  --starting-smiles input.smi --activity activity.txt \
  --found found_sorted.smi --notfound notfound_sorted.smi \
  --predictions notfound.pred --products-sorted -o transformations.html
```

`--products-sorted` trusts that each product file is sorted in increasing atom
count. It retains the first capped candidates per reaction from each file and
compares them to choose the smallest examples overall. Statistics still include
all scored records. Within atom-count ties, the capped candidate set can depend
on the sorted file's tie order.

Several products may share one starting-ID/reaction identifier. Their predictions
can still contribute to the reaction distribution, but cannot safely be assigned
to individual structures. Such examples explicitly show that ambiguity and omit
pair-specific predicted activities. When all corresponding structures are
identical according to LillyMol, their mean predicted activity can be shown.
Scoring records need not follow product-file order.

## Distributions and controls

Each row has two aligned tracks on the same horizontal scale, shared across the
entire report. Experimental values use dots for up to 150 records, and a histogram
for larger groups. Predicted values use a histogram. All histograms use the same
bin boundaries; each track's height is normalized to its own peak. Compare the
printed counts, not bar heights, when assessing evidence volume. Vertical markers
show medians; the dashed line marks zero change.

Large prediction groups are embedded as aggregate histograms rather than lists
of every value. Exact counts, means, medians, quartiles, and ranges use all values.

Controls provide:

- Independent minimum experimental and predicted record counts.
- Both requirements, or either active requirement. A minimum of zero disables
  that requirement in either mode.
- Sorting by mean, median, quartiles, directional fractions, record count, or
  number of distinct starting molecules, using either evidence source.
- Higher-first and lower-first sorting, so no preferred activity direction is
  assumed. Rows lacking the selected metric sort last.
- Transformation-name search.
- Shortlist checkboxes and CSV export. The shortlist persists when filters change;
  export includes all shortlisted rows, including temporarily hidden ones.

`--neutral-width 0.1` treats changes from -0.1 to +0.1 as neither increasing nor
decreasing when calculating directional fractions. It does not change means,
medians, or distributions. The default is zero; choose a width appropriate to the
assay's measurement precision. The neutral width is fixed when generating the
report. Use `--title` to name the report.

Smallest examples illustrate chemistry; they are not necessarily representative
of the largest effect. Neither model agreement nor a small set of experimental
pairs establishes that a transformation will improve a new compound.

## Checks

From `contrib/bin`, run:

```shell
../../run_python.sh -m unittest discover -s test -p medchem_transform_activity_html_test.py
```

Tests cover statistics, example selection and caps, matching and ambiguity,
LillyMol structural identity, input protection, and offline HTML escaping.
If Node.js is available, they also exercise the embedded sorting, filtering,
expansion, shortlist, and CSV controls with a small DOM harness.
