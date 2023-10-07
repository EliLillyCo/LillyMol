# Evidence

When we collect a biological result, we usually end up with a numerical value
associated with each molecule. This data is much more than an individual number
associated with each molecule, there is a story in there. There is evidence.

If we see a particular result for a molecule, and we also see very similar
results for its close neighbours, that is very strong evidence that the data
associated with our starting molecule is good. It is backed up by supporting
evidence.

If on the other hand, we see a particular result for a molecule, but
we see a very different results for its close neighbours, that raises
questions about the validity of the molecule under consideration. It it
an activity cliff? Is it an experimental artifact or a data problem? Is
the structure correct? There is no supporting evidence for the activity
associated with that molecule.

We also frequently encounter the case where a molecule may not have close
neighbours, but moderately close neighbours, and where the neighbours
cannot really provide strong evidence one way or another, for the activity
of the original molecules.

## Evidence
The `evidence` tool reads an activity file and a nearest neighbours file
and, for each molecule, gathers statistics about the molecule, its activity
and the activity of its neighbours. The tool is designed for use in an
interactive structure examination tool like Spotfile, DataWarrior or similar.

In addition to reporting minimal statistics about the nearest neighbours,
the tool calculates a number of pseudo knn 'models' for each target molecule,
and reports those.

## Usage
The usage message is
```
Computes measures of internal consistency for a series of measured values based on how
consistent the values are across the neighbours
Takes a single argument, a TFDataRecord serialized nnbr::NearNeighbours protos, such as
might be produced by the -T option to nn2proto
The following arguments are recognised
 -config <fname>     an EvidenceData::Options textproto file with options
 -A <fname>          descriptor file containing activity values for each molecule
 -smiles <fname>     smiles for each molecule, will be included in the output
 -C                  data is classlfication type (not implemented)
 -diff               for each column generated, insert an extra column with difference from actual
 -v                  verbose output
```

Unlike other tools, this options for this tool are largely driven by the contents
of the textproto configuration file provided by the `-config` option. See the
proto definition to see what options you might like. For testing I have found
this configuration useful.
```
knn: [1, 2, 5, 10]
closest_value: [1, 2, 5, 10]
piecewise_linear {
  min: 0.15
  max: 0.4
}
```

This builds 4 different KNN models, reports the closest activity to the
target within the 1, 2, 5, 10 neighbours, and explores a piecewise linear
weighting function.

The `-A` option is mandatory and is a descriptor file containing the activity data.

If you would like smiles in the output, add the `-smiles` option to provide
a smiles file with the smiles for every identifier.

The tool will ultimately support classification data, but that is not there yet.

If the `-diff` option is specified, for every calculation that the tool does, it will
add a column to the output which contains the (signed) difference between the computed
value and the actual value. This should be useful when sorting in order to identify
outliers.

## Typical Workflow

During testing and development, I used this sequence of commands.
```
gfp_make -STD data.smi > data.gfp
gfp_nearneighbours_single_file -n 10 data.gfp > data.nn
nn2proto -T Tfile -v data.nn
evidence -smiles data.smi -v -A data.activity -config evidence.textproto Tfile
```
