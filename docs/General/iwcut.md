# iwcut

`iwcut` is very much like the standard Linux cut command. It was built
before `cut` adopted some of the features built into `iwcut`. The most
notable remaining difference is that `cut -f 2,1` will write out columns
1 and 2, whereas `iwcut -f 2,1` will write out column 2 followed by
column 1. This makes it very useful for doing things like extracting
a smiles file from something like
```
methane,C,...
```

It also understands column headings, so if you have a tabular file with a
header, it can extract columns by name, while also preserving the column
order requested on the command line. In that case, it also knows that the
first column is special, the identifier column, and that will also be written.

There are a variety of other options specific to the kinds of files we
handle.

One common operation is to fetch the columns that are in one file from
another file.

For example if 'file1.dat' contains
```
id F1 F2 F3
id1 1.1 1.2 1.3
```
and 'file2.dat' contains
```
id F1 F3 F4
id2 2.1 2.3 1.4
```
then asking to extract all the descriptors from 'file1.dat' that are
in 'file2.dat' looks like
```
iwcut -D file1.dat -k . file2.dat
```
which generates
```
id F1 F3
id2 2.1 2.3
```
Descriptor 'F2' was not in 'file2.dat' so with the '-k' option, it was
ignored. Feature 'F4' in 'file2.dat' was not requested.

There are limits to what iwcut can do - it is a fundamentally simple tool.
Using a data analysis package like Julia or Python provides more robust
handling for DataFrame type operations.


