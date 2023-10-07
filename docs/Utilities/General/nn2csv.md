# Near Neighbour Files to CSV

## Background

Most of the gfp_* tools generate output files in a TDT form,
which can then be passed to nplotnn for conversion to smiles or
other forms.

For conversion to smiles forms, it works well, and can easily
perform rudimentary filtering on the neighbour list.

Over time, more and more people have wanted to be able to convert
to tabular form, and I have tried to make nplotnn generate tabular
output. That was a mistake, it really cannot do that and I will
try to remove those attempts.

## nn2csv
This is a very simple tool that is designed to do one thing, and to
do it very well. It converts a nearest neighbour file to csv form.
It hardly does any filtering, and has only a couple of options.
I thought hard about even adding the options I did, perhaps it should
have none.

Like nplotnn, it can process multiple files. Clearly in order to generate
a tabular file, it must scan them all in order to find what is the
maximum number of neighbours, which then controls how many columns
are needed. Molecules that do not have enough neighbours are padded with
'*' characters.

Here is the usage message
```
Converts a .nn file to csv form
 -o <sep>      set token separator (default ,)
 -n <nbrs>     only write <nbrs> neighbours
 -s ...        sort the targets, enter '-s help' for info
 -z            do not write molecules with no neighbours
 -v            verbose output
```
which is refreshingly brief.

with no options, it write all data from the input file. I anticipate
this being the most common use.

Since it was easy, I added the ability to trim the number of columns
generated with the `-n` option. And also since it was easy, I added
sorting - even though such sorting will be trivial once the csv
data is read into `Julia`, `Python`, `R`, `Excel` or other data
analysis tool.

Note that if you ask for more columns than are in the input, it
will dutifully generate the number of columns of output you
requested! The extra columns are empty.

The -z option might save you a subsequent step in filtering, and
it was very easy here.

Same with the sorting option, it was easy to do here, and might
save a step in a slower language.

Most of the other functionality in nplotnn simply does not apply,
the tools are different.

## Pipelines
Note that while this tool does scan all its inputs before generating
any output, it accumulates all the data it needs, and does not
rescan the input. So if needed, it can be used as a pipelined
command. Therefore this command combination is ok
```
gfp_nearneighbours_single_file.sh -v -T 0.5 -n 20 ~/rand.chembl.gfp | nn2csv.sh -
```

