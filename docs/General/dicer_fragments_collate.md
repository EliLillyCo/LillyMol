# dicer_fragments_collate

This tool aggregates multiple `dicer_data::DicerFragment` protos
into an aggregated form. This problem might arise in a situation
where dicer has been run across multiple sets of molecules, each
of which generated its own set of protos, which need to be
aggregated across those multiple results.

For example if one file contained
```
iso: ATT smi: "O[1CH]=O" par: "CHEMBL1213530" nat: 3 n: 100 
```
and another contained
```
iso: ATT smi: "O[1CH]=O" par: "SIGMA28760326" nat: 3 n: 50 

```
the final result would be
```
iso: ATT smi: "O[1CH]=O" par: "CHEMBL1213530" nat: 3 n: 150 
```
where the `n` attribute is the sum of the individual vales,
and the `par` value is that of the first one encountered as
the files are scanned.

## Scalability
When run on large collections, with few restrictions on what fragments
get formed, this tool can consume large amounts of RAM, since the
data is loaded into internal hashes.

## Arguments
The arguments are
```
Aggregates multiple dicer_data::DicerFragment text proto files
 -p <support>        minimum support level (n: value) for inclusion
 -nosmi              each record does not contain a leading non-unique smiles
 -r <n>              report progress every <n> items read
 -minimal            extract only the essential information from the protos
 -tfdata             data is TFDataRecord serialized protos
 -v                  verbose output
```

### -p <support>
Impose a support requirement for output. If the sum of all values is less
than <support> the result will not be written.

### -nosmi
The input is textproto form, and does **not** contain a leading smiles string.

### -r <n>
Report progress every <n> items processed. This can be a long running task and
knowing how it is progressing can be helpful.

### -minimal
The default mode is to store the protos in memory. With this option, only the minimal
amount of data is extracted from the protos. This can help with memory consumption.

### _tfdata
The input data is TFDataRecord serialized protos. The output is always textproto form.
