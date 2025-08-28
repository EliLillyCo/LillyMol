# unique_rows
The standard Linux tool uniq can be used for detecting duplicates in a tabular file. But
the input must be sorted. `unique_rows` does not need to have its input sorted, it maintains
a hash of what has been encountered, so may be more convenient to use.

It may also be more efficient than sorting and using uniq. For example sorting a file
of Chembl smiles takes 1.4 seconds - which seems to have included some multi-threaded
processing. Then running uniq on that takes an additional 0.4 seconds. `unique_rows` does
the whole calculation in about 1.0 seconds, with the advantage of being one command
rather than two, and preserving the order of the input file..

The following options are recognised
```
 -c <col>       only consider column(s) <col>
 -x <col>       consider all column(s) except <col>
 -P <fname>     read a file of identifiers to be considered previously seen - added to the hash before reading the input file.
 -D <file>      file for duplicate records - by default duplicate records are discarded.
 -o             display counts of how many times each identifier encountered
 -O <fname>     write counts and times each identifier found to <fname>
 -s             sort output from both -o and -O. Simlified output too, just token and count
 -z             remove leading 0's from comparisons
 -n             no output, just collect statistics
 -nc            ignore case during comparisons
 -r <number>    number of instances of each unique item to write (default 1)
 -p <prob>      discard duplicate records with probability <prob>
 -trunc <char>  truncate fields at first <char> - useful for truncating decimals
 -tab,-csv      deal with differently formatted files
 -i <...>       input column separator
 -nwz           No Warnings about Zero length text comparisons
 -subset        regexp specifying which part of the column to consider '(..*)...' would be all except last 3 chars
 -h <n>         first <n> records in file
 -j             treat as descriptor file
 -v             verbose output
```
The most commonly used option is the `-c` option, which specifies the column(s)
to consider. For the Chembl case above
```
unique_rows -c 2 -v chembl.smi
```
which (of course) reported no duplicates.

Most of the options are straightforward. The -O option can be useful, generating a
summary of tokens and how often they are encountered. This can be useful in some
situations.
