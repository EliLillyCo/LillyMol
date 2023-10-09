# fetch_smiles_quick

Extract selected rows from a file. Performs an inner join between columns
in two different files. The most common usage is to start with two files

1. File containing identifiers of interest
2. File containing a superset of the identifiers of interest.

```
fetch_smiles_quick idfile file.smi > idfile.smi
```
where the identifiers in `idfile` are assumed to be in column 1 and
the identifiers in the smiles file, `file.smi` are assumed to be
in column 2.

Note that the order in which the output is written is according
to position in `file.smi`. If it is important to preserve the
order specified in `idfile`, use `fetch_smiles` instead.

# HOWTO
The following options are recognised.
```
Fetches records from one file based on identifiers in one or more other file(s)
 -c <col>       identifier column in identifier file
 -C <col>       identifier column in smiles file
 -C RX=<rx>     identifier is whichever column(s) match <rx>
 -d             ignore duplicate identifiers in identifier file
 -q             quietly ignore duplicate identifiers in identifier file
 -a             write all instances of identifiers in smiles file
                by default, only the first is written
 -X <fname>     write smiles records not in <identifier_file> to <fname>
 -Y <fname>     write identifiers not in <smiles_file> to <fname>
 -w             write the -Y file as a smiles file (swap columns)
 -k             suppress addition of info from identifier file
 -n <string>    string to insert between record and info from identifier file
 -x             invert behaviour, selection becomes deselection
 -z             strip leading zero's from identifiers
 -j             identifier file is descriptor file, skip header record
 -b             stop processing identifier file on error, but continue processing (dangerous)
 -i <char>      column separator in identifier file
 -I <char>      column separator in smiles file
 -S <stem>      first files are identifier files, last is haystack. Create many subsets
 -u <suffix>    suffix for -S files created
 -g <ndx>       start number for files created (-g 1 for dopattern)
```

# Options

## -c_\<col\>
Identifier column in identifier file. Default is column 1.

## -C_\<col\>
Identifier column in smiles file. Default is column 2.

## -C RX=\<rx\>
Identifier is whichever column(s) match \<rx\>.

## -d
Ignore duplicate identifiers in identifier file. By default duplicate
identifiers in the identifier file is a fatal error.

## -q
Quietly ignore duplicate identifiers in identifier file.

## -a
Write all instances of identifiers in smiles file - by default,
only the first is written.

## -X \<fname\>
Write smiles records not in \<identifier file\> to \<fname\>.

## -Y \<fname\>
Write identifiers not in <<smiles file\> to \<fname\>.

## -w
Write the -Y file as a smiles file (swap columns).

## -k
Suppress addition of info from identifier file.

## -n \<string\>
String to insert between record and info from identifier file.

## -x
Invert behaviour, selection becomes deselection.

## -z
Strip leading zero's from identifiers.

## -j
Identifier file is descriptor file, skip header record.

## -b
Stop processing identifier file on error, but continue processing (dangerous).

## -i \<char\>
Column separator in identifier file. Note that things like `-i tab` are
recognised, and converted to their non printing form. Enter `-i help` for
a list of the directives that are recognised.

## -I \<char\>
Column separator in smiles file.

## -S \<stem\>
First files are identifier files, last is haystack. Create many subsets.

## -u \<suffix\>
Suffix for -S files created.

## -g \<ndx\>
Start number for files created (-g 1 for dopattern).
