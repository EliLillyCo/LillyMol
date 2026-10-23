# mispredicted correctness review

The implementation is `src/Utilities/General/mispredicted.cc`.
This records the review findings so they can be addressed in subsequent changes.

## Addressed findings

- Bound worst-count selection by the number of predicted molecules.
- Bound balanced selection, including single-molecule and one-sided data.
- Scan threshold selection safely when all or no molecules qualify.
- Activate valid filename regexes and reject invalid ones.
- Reject invalid numeric predictions in normal mode; skip invalid prediction
  records and continue in fault-tolerant mode, retaining valid observations.
- Keep the first experimental value with `-h`, including classification inputs.
- Reject summary/textproto paths that alias inputs or each other before opening
  outputs. This includes inputs from file lists and filename regexes.
- Return failure when there are no valid predictions or primary output fails.
- Reverse continuous output with `-r` after filtering, including summary and
  textproto output, without changing selected molecules.
- Break remaining ranking ties by ID, including classification output.

Balancing uses accumulated error magnitude in each direction while both signs
remain available, and writes at most the requested count. Continuous output and
`-w`/`-W` now use mean absolute error by default. Explicit `-s` selects mean
signed error (observed minus predicted); filters use its magnitude. `-b` requires
`-w` and implies signed errors to preserve directional balancing.

Single-prediction summary and textproto records report zero standard deviation
without calling variance. Multiple predictions use the sample standard deviation
of the selected error metric.

## Remaining findings

- Classification accepts `-w`, `-W`, `-Y`, and `-T` without implementing the
  corresponding selection or output. `-S -d` adds a SMILES header but no SMILES
  values in classification rows.
- The classification textproto writer does not populate the class label.
- `-D` is advertised and accepted but never processed. Support-vector distances
  are stored with `-k`, but their intended output and binning are unused.
- Equal class prevalence makes `-X` class translation depend on hash iteration
  order; explicit tie handling is needed.
- Duplicate prediction records are counted separately; independent-run validation
  remains a possible follow-up.
- Output error propagation in summary/textproto helpers remains incomplete.
- Continuous input does not translate tabs, unlike classification input.
- Entity and per-class accumulator ownership uses leaking raw pointers.

Fault tolerance should retain valid records and ignore invalid records, rather
than discard entire files. This is the intended behavior, not a requirement for
transactional file processing.

## Minimum prediction coverage

`-n <count>` requires a positive number of valid prediction records per identifier
(default 1). Regression and classification apply this restriction before ranking
or selection. All implemented output paths use the eligible identifiers. A filter
that leaves no eligible identifiers succeeds with empty data output (headers are
retained where applicable). An input with no valid predictions remains an error.
Repeated records within one file count separately; the minimum is not a count of
independent files or models. Invalid records skipped with `-F` do not count.
