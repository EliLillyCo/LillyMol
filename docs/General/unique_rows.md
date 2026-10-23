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
 -column <spec> independent per-column limit, e.g. '1,n=4,subset=(.*)..$'; repeatable
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

## Independent limits for multiple columns

Use repeated `-column` specifications to limit each column independently:

```shell
unique_rows -column '1,n=1' -column '2,n=88,subset=(.*)..$' input.txt
```

A row is retained only if every configured column is below its limit. Counters
are incremented only for retained rows: rejection by one column does not consume
another column's allowance. Input order determines which rows are retained.
Counters persist across input files.

Each specification starts with a one-based column number; negative numbers count
from the end (`-1` is the last column). Optional `n=` sets a positive limit,
defaulting to 1. Optional `subset=` extracts the single captured group from a
regexp that must match the entire column value. Omit it to compare the whole
value. `subset=` must come last, allowing commas within the regexp. Invalid
specifications, missing columns, and regexp non-matches are errors.

For example, given:

```text
a x
a y
b y
c y
d y
```

`-column '1' -column '2,n=2'` retains `a x`, `b y`, and `c y`.
The rejected `a y` does not consume a `y` allowance.

This is separate from `-c 1 -c 2`, which retains its existing behaviour of
combining columns into one uniqueness key. `-column` cannot be combined with
`-c`, `-x`, `-R`, `-d`, `-r`, `-subset`, `-o`, `-O`, `-P`, `-p`, or `-all`.
Repeated specifications for the same column number are rejected.

Input separators, quoted fields, `-h`, duplicate output (`-D`), and suppression
of normal output (`-n`) remain available. `-z`, `-nc`, and `-trunc` operate on
each selected value before its optional subset regexp. As in legacy mode, a
value starting with the truncation character is left intact with a warning.
Empty values are counted
as ordinary values in this mode. Use shell quoting around specifications
containing regexp metacharacters.

## Preparing comparison values

Comparison operations change the keys used for uniqueness, not the records
written to output. They run in this order:

1. `-z` removes leading zeros from each selected token.
2. In legacy mode, selected tokens are joined with spaces into one comparison
   key. With `-column`, each selected column remains an independent key.
3. `-nc` converts comparison text to lowercase.
4. `-trunc <char>` removes the first occurrence of that character and everything
   after it. If the text starts with the character, a warning is issued and the
   text is left intact.
5. The optional subset regexp extracts its captured group. Legacy `-subset`
   applies to the combined key; `subset=` applies to its individual column.

The regexp must match the entire comparison text after the preceding operations.
For `-column`, exactly one capture group is required. For example:

```shell
unique_rows -column '1,subset=(.*)..$' -nc input.txt
```

The values `ABC01` and `abc02` both produce the comparison key `abc`, so only the
first row is retained with the default limit of one. The original capitalisation
and suffix remain in the output row.

To ignore decimal suffixes independently in two columns:

```shell
unique_rows -column '1,n=2' -column '2,n=3' -trunc . input.txt
```

Here `12.5` and `12.8` compare as `12`. A value such as `.5` is left intact with
a warning, preserving the existing truncation convention.
