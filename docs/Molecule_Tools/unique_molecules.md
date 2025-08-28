# Unique Molecules

## Background
A common task in Cheminformatics is the identification of unique or
duplicate structures. The unique smiles is commonly used for this.

`unique_molecules` reads a stream of molecules, and if any have been
seen before, that molecule is discarded.

## HOWTO
The usage message is
```
  -l             strip to largest fragment
  -a             compare as tautomers - skeleton and Hcount
  -c             exclude chiral info - optical isomers will be duplicates
  -z             exclude cis/trans bonding information
  -I             ignore isotopic labels
  -y             all non-zero isotopic values considered equivalent
  -f             function as filter (TDT input)
  -p <fname>     specify previously collected molecules
  -G <tag>       identifier tag when working as a filter
  -s <size>      specify primary hash size (default 1000)
  -S <name>      specify output file name stem
  -D <name>      write duplicate structures to <name>
  -R <rxn>       perform reaction(s) on molecules before comparing
  -T             discard molecular changes after comparison
  -r <number>    report progress every <number> molecules
  -e             report all molecules together with counts
  -d             gather data as DicerFragment protos, -U file is more informative
  -U <fname>     write molecules and counts to <fname>, add '-U csv' for csv
  -d             use '-U smiles' to write smiles + textproto
  -j             items are the same only if both structure and name match
  -t E1=E2       element transformations, enter '-t help' for details
  -i <type>      specify input type
  -o <type>      specify output type(s)
  -A <qualifier> Aromaticity, enter "-A help" for options
  -K ...         standard smiles options, enter '-K help' for info
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -v             verbose output

Typical usage might be

unique_molecules -l -c -z -I -S /tmp/unique -v /home/ian/rand50k.smi

To generate DicerFragment protos in the -U file.
unique_molecules -l -c -z -I -U Ufile.textproto -d -U smiles -v /home/ian/rand50k.smi
```

There are a great many ways by which two molecules can be considered identical.
`unique_molecules` provides a wide variaty of structural modifiers that can
be applied to molecules before their unique smiles are generated.

* -l strip to the largest fragment
* -c discard chirality
* -z discard cis/trans bonding information
* -I discard isotopic labels.
* -y all non zero isotopic labels compared as equivalent
* -R \<rxn\> apply a transformation reaction to the molecules.
* -t E1=E2 transform all elements E1 to E2 (suggest -T I=Cl -T Br=Cl)

Once these transformations are applied, the unique smiles is generate
and compared to what has been encountered previously. If this is
the first instance of the smiles, that molecule is written to the output
strea, otherwise it is classified as a duplicate, and if specified,
written to the -D file for duplicates.

By default, the changed molecule will be written, but the original form
can be written if the -T option is used. That should probably be the
default.

If there is a previous set of molecules (`-p`) that should be considered,
which will establish the unique smiles hash, without writing anything.

By default, smiles are accumulated in a C++ set, which has only keys. If
you wish to get a summary of smiles and counts, use the -e option,
which will use a map for the unique smiles, and counts can be accumulated.

# Options
Some of the other options modify how this works.
## -p \<fname\>
Specify a file of previously selected molecules. The unique smiles hash
is first populated with these molecles, and then the input is compared.

## -G \<tag\>
Obscure.  Identifier tag when working as a TDT filter - GFP files.

## -s \<size\>
Ignore. Specify primary hash size (default 1000). Usually not needed these
days.

### -S \<name\>
Specify output file name stem.

### -D \<name\>
Write duplicate structures to \<name\>. By default, duplicate structures
are just discarded.

### -T
Discard molecular changes after comparison. This is important. For example
if chirality is discarded during unique smiles comparison, by default,
there will be no chirality in the ouput. With the -T option, the original
molecule is retained. This should be the default.

### -r \<number\>
Report progress every \<number\> molecules processed. Or just look
at what is getting written.

### -e
It can be helpful to get a summary of unique smiles and counts 
written to the log.

### -U \<fname\>
Formalise what the `-e` option does, writing the data to a
separate file - rather than mixed in with stderr.

### -j
Only consider structures the same if their names also match.

# Large Datasets
Some thoughts on handling larger datasetsdatasets.

A tool like unique molecules must necessarily scale with the
number of molecules in the input, and the memory requirments
will scale with the number of distinct structures encountered.

Up to several tens of millions of molecules, this will be OK.

For example, running unique_molecules on 34M molecules generated
by making changes to Chembl, runs in about 29 minutes, while consuming
as much as 7GB of RAM. This may be just fine.

The process can be sped up in a number of ways.

The most straightforward is to split the dataset into different
streams, that are guaranteed to be different. For example if the
input is split into files, each of which has molecules with a given
number of heavy atoms, those files can be processed separately. Clearly
a molecule with 10 atoms need not be compared with a molecule with 11 in 
order to find duplicates.

Obviously these separate files can be processed in parallel, and the memory requirements
for each job will be modest.

Taking this concept one step further, perhaps the file could be
divided by molecular formula, with files containing only molecules
with the same molecular formula. Again, files can be processed independently
and memory requirements will be lower.

# A Failed Idea.
Rather than creating separate files, how about we just sort the file
by molecular formula and write something to do uniqueness, but aware
of the molecular formula from one line to the next - and never compare any two molecules with
different formulae.

Given the set of 34M molecules, we can run
```
fileconv -g all -p noname -p MF -S with_formula file.smi
```
which appends the molecular formula to the names in `with_formula`.
```
NC1=CC=C(N(CCF)CCF)C=C1 CHEMBL148945 C10N2F2H14
```
This can be sorted with standard unix `sort`. This is done in about 25 seconds,
but consumes 7.6 GB of RAM.
```
sort -k 3,3 --parallel 8 ...
```
`sort` is very sophisticated, and has options for specifying the buffer
size to use. For larger datasets, this could be used for limiting the
RAM used, and using temporary files.

Once the sorted file is available, the tool `unique_molecules_sorted` can
be used to look for duplicates. If the molecular formula is different from
one row to the next, smiles interpretation is not needed. But if there
are consecutive lines with the same formula, smiles interpretation is
needed. If that is done, the file can be processed in just under 12 minutes.
This could of course be overlapped
```
sort -k 3,3 -parallel 8 file.smi | unique_molecules_sorted -C 3 -v -
```
which would avoid the need for a large temporary file.

Overall, this is not worth it. We can process the original file in 29
minutes using `unique_molecules`. The alternative is a three step process
1. 13 minutes - use fileconv to append the molecular formula
2. 25 seconds - sort the resulting file
3. 12 minutes - unique molecules using the molecular formula as a guard.

For a total run time of about 26 minutes. Not much different from the
29 minutes used by unique_molecules. And the highest RAM requirement
is about the same - this time with 'sort'.

The best way to do approach this problem is likely to create separate
files, split by formula, and process those separately and in parallel.
There is no existing LillyMol tool to do this, but it would be easy to 
do.

## Inchi
Note that if LillyMol has been built to include Inchi functionality,
unique_molecules will respond to the -C option by doing all comparisons
via the Inchi key. This can be quite expensive, but may yield important
equivalence relationships that would not otherwise happen.
