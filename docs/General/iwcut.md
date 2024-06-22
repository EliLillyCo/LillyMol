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
first column is 'special'.

There are a variety of other options specific to the kinds of files we
handle.
