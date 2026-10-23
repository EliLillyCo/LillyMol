# Integer Descriptors to Fingerprints

This directory contains example files for `descriptors_to_fingerprint_integer`.

This example demonstrates how integer descriptor vectors can be converted to GFP fingerprints while preserving similarity relationships. The included data set contains a prototype vector, several closely related variants, and a collection of unrelated vectors.

Generating the example data requires Julia.
The generated CSV files are included in the repository and regeneration
is not normally necessary.

For details on command line options and fingerprint formats, see
[Main Documentation](../descriptors_to_fingerprint_integer.md)

The .csv files contain randomly generated data from ChatGPT. The julia code
to generate the files is `vector_similarity.jl`. That file generates a
prototype vector, then 10 minor variants of that vector. The rest of
the dataset are randomly generated. The variants of `prototype` are labelled
"active" and the randomly generated as "inactive".

Running

```bash
julia vector_similarity.jl
```

generates:

```
  prototype.csv
  synthetic.csv
  prototype_similarity.csv
```

The file `prototype_similarity.csv` shows Tanimoto similarity computations
performed in Julia 
can be used to verify that similarity values computed by LillyMol agree with
independently calculated reference values.

Convert those files to fingerprint form with arbitrary smiles. Use the `-W`
option to create a separate fingerprint file for each .csv input file and
place those files in the current directory, `-W .`
```bash
descriptors_to_fingerprint_integer -W . -s C -i , prototype.csv synthetic.csv 
```
generates `prototype.gfp` and `synthetic.gfp`.

Use gfp_lnearneighbours to find the nearest neighbours of `prototype.gfp` in
`synthetic.gfp`
```bash
gfp_lnearneighbours -p prototype.gfp -n 20 synthetic.gfp > prototype.nn
```
The first ten neighbours should be the records labelled "active".
This demonstrates that the sparse fingerprint representation preserves
the similarity relationships present in the original descriptor vectors.

```bash
grep '^PCN' prototype.nn
```
should show 'active0, active1, active2...`, not necessarily in that order,
as the closest neighbours.
```

That .nn file can be processed with [nplotnn](/docs/GFP/nplotnn.md).
