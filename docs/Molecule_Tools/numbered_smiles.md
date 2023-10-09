# numbered_smiles

## Purpose
`numbered_smiles` is a tool intended for use by expert. It consumes a structure
file and writes output that will have either isotopic labels, or atom map numbers
corresponding to whatever kind of numbering is specified on the command line. This
can be very useful for investigating complex structural issues.

In its most simple form, given ethane, `CC`
```
numbered_smiles ethane.smi
```
will generate a file `ethane_num.smi` that contains `C[1CH3] ethane`, where atom
0 has been assigned isotope 0, no isotope, and atom 1 has isotope 1.

`numbered_smiles` is used for the `-N` option to vf.

## Options
The help message is
```
  -n atom        number by atom number (the default)
  -n symm        number by symmetry class
  -n canon       number by canonical order
  -n ring        number by ring system
  -n <number>    number by atom number starting with <number>
  -e <step>      number step (default 1)
  -P <atype>     number by atom typing (no -n option needed)
  -I             overwrite any existing isotopic information
  -r             sort molecule by numbering
  -p             change to graph form (all bonds become single etc...)
  -H             make all implicit Hydrogens explicit
  -b ...         write bond symmetry information to stderr
  -q <query>     only label atoms matched by the query/queries
  -s <smarts>    only label atoms matched by the smarts
  -m             apply atom map numbers instead of isotopic labels
  -S <stem>      specify output file name stem
  -K ...         standard smiles control options
  -i <type>      specify input type
  -A <qualifier> Aromaticity, enter "-A help" for options
  -v             verbose output
```


