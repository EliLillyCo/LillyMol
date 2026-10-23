# Unique Molecules

## Quick Start

Check for duplicates (ignoring chirality) in a set of molecules.

```bash
unique_molecules -c input.smi
```

Remove duplicate molecules while ignoring chirality, fragments and isotopes.

```bash
unique_molecules -v -l -c -I -S unique input.smi
```
duplicate molecules are discarded and unique molecules written to 'unique.smi'.

## Background
A common task in Cheminformatics is the identification of unique or
duplicate structures. A hash of unique smiles is commonly used for this.

`unique_molecules` reads a stream of molecules, performs any requested
transformations on the molecule and forms a canonical smiles of the
possibly transformed molecule.  If that smiles has been encountered
before, that molecule is discarded, leaving only unique molecules in
the output stream. The new unique smiles is added to the internal hash.

## Limitations

`unique_molecules` keeps all unique structures in memory.
For extremely large datasets, use [unique_molecules_parallel](unique_molecules_parallel.md)
or partition the dataset first.

## HOWTO
The usage message is

## Structural Equivalence Options
```
 -l             strip to largest fragment.
 -a             compare as graph/tautomers - concatenate molecular skeleton and formula.
 -c             exclude chiral info - optical isomers will be duplicates.
 -z             exclude cis/trans bonding information - true by default.
 -I             ignore isotopic labels.
 -y             all non-zero isotopic values considered equivalent.
 -R <rxn>       perform reaction(s) on molecules before comparing.
 -t E1=E2       element transformations, enter '-t help' for details.
 -j             items are the same only if both structure and name match.
 ```
 ## Output Options
 ```
 -S <fname>     specify output file name stem - the unique molecules.
 -D <fname>     write duplicate structures to <fname>.smi.
 -U <fname>     write molecules and counts to <fname>.
 -U csv         write the -U file a csv form.
 -U textproto   write the -U file as dicer_data::DicerFragment textproto
 -U smiles      when writing textprotos, write smiles + textproto
 ```
 ## Other Options
 ```
 -p <fname>     specify previously collected molecules.
 -T             discard molecular changes after comparison - writes initial molecule.
 -r <number>    report progress every <number> molecules.
 ```
 ## Common LillyMol options
 ```
 -i <type>      specify input type. All LillyMol input formats supported, including stdin.
 -o <type>      specify output type(s).
 -E ...         standard element options, enter '-E help' for info.
 -A ...         standard aromaticity options, enter '-A help' for info.
 -K ...         standard smiles options, enter '-K help' for info.
 -g ...         chemical standardisation options, enter '-g help' for info.
 -v             verbose output.
```

There are a great many ways by which two molecules can be considered identical.
`unique_molecules` provides a wide variety of structural modifiers that can
be applied to molecules before their unique smiles are generated.

* -l strip to the largest fragment.
* -c discard chirality.
* -z discard cis/trans bonding information - true by default.
* -I discard isotopic labels.
* -y all non zero isotopic labels compared as equivalent
* -R \<rxn\> apply a transformation reaction to the molecules.
* -t E1=E2 transform all elements E1 to E2 (suggest -t I=Cl -t Br=Cl)

Once these transformations are applied, the unique smiles is generated
and compared to what has been encountered previously. If this is
the first instance of the smiles, that molecule is written to the output
stream, otherwise it is classified as a duplicate, and if specified,
written to the -D file for duplicates.

By default, the changed molecule will be written, but the original form
can be written if the -T option is used. That should probably be the
default. This only applies if the -l, -t or -R options are used.

It is important to note that unique_molecules operates by discarding
duplicates. If the -c (ignore chirality) option is present, and there
are chiral duplicates in the input, the *first* variant is what will
be written. If you need to identify the actual duplicate molecules
use [common_names](common_names.md).

If there is a previous set of molecules (`-p`) that should be considered,
which will establish the unique smiles hash, without writing anything.

By default, smiles are accumulated in a C++ set, which has only keys. If
you wish to get a summary of smiles and counts, use the -U option,
which will use a map for the unique smiles, and counts can be accumulated.

# Options
Some of the other options modify how this works.
## -p \<fname\>
Specify a file of previously selected molecules. The unique smiles hash
is first populated with these molecules, and then the input is compared.
All structure modification options are applied.

### -S \<name\>
Specify output file name stem.

### -j
Only consider structures the same if their names also match.

### -T
By default, any transformations applied for duplicate comparison are retained in the output molecule.

Use `-T` to restore the original molecule before writing output.
This is often desirable when using `-c`, `-I`, `-R`, or `-t`.

### -D \<name\>
Write duplicate structures to \<name\>. By default, duplicate structures
are just discarded. The default output of the -D option file is
just a smiles file with the smiles and id of all duplicates encountered.

### -U ...
At the end of processing write a file containing summary information
on the internal hash structures. By default, this will be a list of
unique smiles and number of times it has been found, '-U fname'.
```
NS(=O)(=O)OC1CCCC1 1
N#Cc1sc(N(=O)=O)cc1 1
O=C(OC)C1C(CC2N(C)C1CC2)c1ccc(Br)cc1 2
Clc1cc(c2ccoc2)c2N3CCCC3NS(=O)(=O)c2c1 2
```
That file can be written as csv form via '-U fname -U csv'.

It can be written as a dicer_data::DicerFragment textproto form by
'-U fname -U textproto' which might look like
```
smi: "N[C@@H]1NCC(C)C([C@@H]1O)O" par: "CHEMBL85218" nat: 10 n: 2 
smi: "Brc1cc(C)c[n+](C)c1Cl" par: "CHEMBL44713" nat: 10 n: 2 
smi: "CN(Cc1ccccc1)C" par: "CHEMBL45591" nat: 10 n: 2 
smi: "O=C(C)CCCN(=O)(C)C" par: "CHEMBL25004" nat: 10 n: 2 
smi: "N=C(N)Nc1[n]cc[n]c1" par: "CHEMBL54990" nat: 10 n: 2 
```
### -r \<number\>
Report progress every \<number\> molecules processed. Or just look
at what is getting written.

### -g ...
Apply chemical standardisation. This is likely very important if you
are dealing with disparate sources of molecules. For example if
one data source writes Nitro groups as 'O=[N+]-[O-]' and the
other writes them as 'O=N=O', the unique smiles of the starting
molecules will be different. If '-g all' is specified, they become
the same. See [chemical_standardisation](/docs/Molecule_Lib/chemical_standardisation.md)

# Algorithm.
Internally the tool uses an array of hash_set's to keep track of
the unique smiles encountered by atom count. The reason for segregating
by atom count comes from the observation that performance seemed
to degrade as the hash became larger. So segregating the dataset
by atom count provided a means of mitigating that problem. It is
not clear whether today's hash implementations suffer from that
issue, but the current implementation should not be deleterious.

For discussion of large-scale parallel processing, see
[unique_molecules_parallel](unique_molecules_parallel.md).

## Example Invocations.
In each case, add '-S output' in order to write the unique molecules
to 'output.smi'. Also, adding chemical standardisation '-g all' is
usually a good idea.


| Task | Command |
|---|---|
| Ignore chirality | `unique_molecules -c input.smi` |
| Ignore isotopes | `unique_molecules -I input.smi` |
| Keep duplicates | `unique_molecules -D dups.smi input.smi` |
| Use prior molecule set | `unique_molecules  -c -p existing.smi input.smi` |
| Match name and structure | `unique_molecules -j input.smi` |

Note too that like most LillyMol tools, unique_molecules usually works
silently. Add the -v option to get diagnostic information about
processing.

## Inchi
Note that if LillyMol has been built to include Inchi functionality,
unique_molecules will respond to the -C option by doing all comparisons
via the Inchi key. This can be quite expensive, but may yield important
equivalence relationships that would not otherwise happen.

