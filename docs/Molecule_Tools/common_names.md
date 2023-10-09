# common_names

## Purpose
There are two utilities that deal with duplicate molecules, `common_names`
and `unique_molecules`.

If all you care about is eliminating duplicates as quickly and efficiently
as possible, use `unique_molecules`.

The purpose of `common_names` is to allow identification of which 
molecules are duplicates.

## HOWTO
The help message for `common_names` is
```
Gather all the names duplicate structures together
  -a             compare graph forms - add 2nd -a option to include H count
  -c             exclude chirality information
  -x             exclude directional bonds
  -l             strip to largest fragment
  -I             remove isotopes before storing
  -T ...         standard element transformation options, enter '-T help'
  -D <separator> separator for when storing duplicate entries
  -f             single pass operation, smiles output only
  -y             write first name and count of smiles only
  -s <size>      maximum number of molecules to process
  -r <number>    report progress every <number> molecules processed
  -i <type>      specify input file type
  -S <name>      specify name for output
  -o <type>      specify output file type
  -E ...         standard element options
  -A <qualifier> Aromaticity, enter "-A help" for options
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -v             verbose output
```
## Background
Molecules can be compared for equality in different ways. Several
command line options govern what components of the input molecules
go into the comparison.

## Options
### -a
Convert molecules to graph form, remove chirality, isotopes and
charges, all bonds become single bonds, and compare the unique
smiles of the resulting molecules. In this way, benzene and cyclohexane
are identical. Add a second `-a` option and the hydrogen count of
the starting molecule is included with the unique smiles, so 
benzene and cyclohexane would no longer be the same.

The molecular graph is useful for doing tautomer comparisons, but
may in some cases be too loose.

### -c
Remove chirality before comparing. Be cognisant of the accuracy
of chirality information you may have.

### -x
Remove cis-trans bonds before comparing. This is highly recommended
since this information is frequently not accurate.

### -l
Strip counterions before comparing.

### -I
Remove any isotopic information before comparing.

### -D \<separator\>
`common_names` generates an output file that has molecules found to
be duplicates grouped together. So if the input contained
```
C methane
C CH4
```
the output file would contain
```
C methane:CH4
```
where the `:` character is set by the -D option.

### -f
If you are just writing smiles, this can be a faster, but memory
requirements may be larger. Definitely worth starting here, and if 
it works, great.

### -s \<size\>
Sets an initial size for certain internal data structures. It is not
used when using the `-f` option, but when processing multiple input
files, it is required. Set to the number of molecules in your
input file(s).

### -r \<number\>
Report progress every \<number\> molecules processed. This can be
helpful when processing large datasets.

### -y
Only works in the absence of the `-f` option. Rather than concatenate
all names found, just write the first one and a count. So, our
methane example would be written as
```
C methane 2
```
### -D \<fname\>
As duplicate structures are identified, write them to \<fname\>.

### -R \<rxn\>
Apply one or more reactions to molecules before doing the comparison.
Not sure if this is a good idea or not, but today only old style
reaction files are supported.

### -T
Frequently it is desirable to consider the heavy halogens to be
equivalent. Use the -T option to enable this.
```
-T 'I=Cl' -T 'Br=Cl'
```
enables such a transformation.

### -S \<stem\>
Specify a file name stem for the output.

The other options are common across different tools.
