# rotatable_bonds
This tool computes and may optionally filter molecules according to
the rotatable bond count.

Invoked with no options, a typical output might look like
```
C(=O)NCCCC CHEMBL45466  ROTBOND = 3
```
where the number of rotatable bonds is appended to the smiles and name.

To filter to a maximum number of rotatable bonds, use the `-M` option
```
rotatable_bonds -v -M 10 file.smi > okrotbond.smi
```

Whether or not a bond is considered rotatable can be a complex decision.
This tool has considerable flexibility in determining what is considered
rotatable. That said, the defaults seem to work well for most situations.
It can also use the default rotatable bond computation that is shared across
several LillyMol tools.

## Options
The following options are recognised
```
  -m <number>    discard molecules with <number> or fewer rotatable bonds
  -M <number>    discard molecules with <number> or more  rotatable bonds
  -w <natoms>    write molecules with <natoms> or fewer atoms regardless
  -f <number>    longest allowable consecutive flexible bonds
  -f <maxc=nn>   maximum connectivity along a long flexible chain
  -f <breaku>    stop any flexible chain growth at any unsaturated atom
  -q <query>     query specification (alternate to -s)
  -s <smarts>    bonds across first two matched atoms are NOT rotatable
  -k             bond is unmatched atom attached to first matched atom
                 can use '-s [CD4](-F)(-F)(-F) -k' rather than
                 can use '-s *-[CD4](-F)(-F)(-F)' which is much slower
  -r <number>    place isotope <number> at the ends of rotatable bonds
  -R <size>      single bonds in rings of size <size> are rotatable
  -a             append rotatable bond count to name
  -c             append longest distance between rings to name
  -C <bonds>     discard molecules with longest distance between rings > <bonds>
  -w <natoms>    write molecules with <natoms> or fewer atoms regardless
  -S <string>    create output files with name stem <string>
  -B <string>    write discarded molecules to <string>
  -V <fname>     list of rotatable bonds written to <fname>
  -F ...         compute rotatable bonds between features - enter '-F help' for info
  -Q ...         use quick and dirty calculation, enter '-Q help' for info
  -e             compute rotbonds on largest fragment only
  -l             reduce to largest fragment (discard counterions)
  -E ...         standard element options
  -i <type>      specify input file type. Enter '-i help' for details
  -o <type>      specify output file type(s)
  -A <qualifier> Aromaticity, enter "-A help" for options
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -v             verbose output
```

The most common operation is filtering a list of molecules to limit
the number of rotatable bonds. The `-m` (lower limit) and `-M` (upper
limit) options control this.

In some circumstances it may be desirable to write small molecules regardless
of rotatable bonds, and the `-w` option does this.

The `-f` option allows finer control of when a molecule is rejected, this
time restricting scrutiny to just specific long chains - rather than an
aggregation of rotatable bonds across regions of the molecule.

The `-q` and `-s` options allow specificatin of bonds that are considered
not rotatable. The query must contain at least two atoms, and the bond
between the first two matched atoms will be considered non rotatable.

The `-r` is a useful diagnostic that allows you to see which bonds it
has considered to be rotatable. Unfortunately it is quite likely you will
encounter cases where a non-rotatable bond will have isotopes at each end
since the adjoining bonds are rotatable. There are ideas for fixing this
if it ever becomes important.

By default, output is to stdout, but that can be changed with the `-S` option,
and if discarded molecules are needed, those can be captured via the `-B` option.

To use the default calculation, used across LillyMol, use the `-Q` option,
and I would recommend `-Q better` to enable a more expensive, but more
realistic computation.

## Between Features
The -F option allows computation of rotatable bonds between features.
This requires specification of two sets of queries, called 'q1' and 'q2'.

For example to compute the number of rotatable bonds between Fluorine
atoms and '[OH]' atoms, that might look like
```
rotatable_bonds -v -F write=/tmp/rotb -F q1:SMARTS:F -F 'q2:SMARTS:[OH]' file.smi
```

Output (/tmp/rotb) might look like
```
C1(=C(C=CC(=C1O)O)[C@H](O)CN)F CHEMBL40796 0
N1(C(=O)N=C(N)C(=C1)F)[C@H]1O[C@@H](CO)OC1 CHEMBL148993 2
C1C(O)(C(F)(F)F)N(N=C1C(C)(C)C)C(=O)N CHEMBL455741 0
C1=CC(=C2NC(=O)C(=C(C(=O)OCC)C2=C1)O)F CHEMBL2029354 0
FC1=CC=C2N(C(C)CCC2=C1)C(=O)CCC(=O)O CHEMBL1488447 1
COC1=C2C(=CC(=C1)N)SC(N)=N2 CHEMBL587076 .
```
The queries are specified the same way that the `-q` option to
`tsubstructure` is handled. Any number of 'q1:' and 'q2:' 
directives can be given. The first query that matches
the molecule will be used.

The output file for this functionality must be specified directly.
All molecules in the input are written to this file. Alternatively
we could send the output to the stream for passing or failing
molecules instead, avoiding the need to specify this here. Let me
know what you would prefer. It might actually be easier to do that...

Molecules that do not match both of the queries will have a '.'
as their value. This file should be easy to sort and filter as
needed. There is no header record, although there could be.

## Other functionality
The `-c` and `-C` options do not really belong in this tool,
since they deal with distances between rings, rather than anything
related to rotatable bonds. To filter molecules to a given distance
between rings, try
```
tsubstructure -s '[R]...{>9;[R0]}[R]' file.smi
```
There must be more than 9 atoms between the two rings, and
none of the atoms on the shortest path between the two ring
atoms are in a ring.
