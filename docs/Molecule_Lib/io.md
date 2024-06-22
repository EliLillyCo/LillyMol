# Molecule I/O

## Concepts
Most LillyMol executables read chemical structures, perform some operation,
and write results in some form. Almost all tools use the `-i` option
to control aspects of the input.

Because almost all tools share the same input scheme, the following applies
to almost all LillyMol tools. New tools are strongly encouraged to use
the same scheme. The file `Molecule_Tools/skeleton.cc` is a prototypical
tool that should be readily adaptable to new tasks.

## Philosophy
### Errors
LillyMol reads molecules until it encounters an error - a molecule
that it cannot read. At that point, the tool will stop with a non zero
error code. Fix your input and try again.

Sometimes however it might be safe to skip over known errors, so there
is a means of telling the tool to ignore a certain number of otherwise
fatal connection table errors. These erroneous molecules are never
given to the tool doing the work - since instantiation of a Molecule
has failed.

### Problems
Many non-LillyMol tools accept unusual/ambiguous input and attempt to convert it
to what the tool thinks you probably wanted. This is completely opposite
the LillyMol strategy which is to assume that the input is exactly
as specified, and to reproduce that as faithfully as possible. There
are mechanisms, which are never default behaviour, to make certain
changes. LillyMol assumes that your input is correct and not to be messed with.

## Usage
Entering `-i help` to almost any LillyMol tool will produce a usage
message
```
 -i sdf                  SDF input
 -i smi                  smiles input
 -i tdt                  TDT input
 -i mdl                  MDL format (generally, use 'sdf' instead)
 -i textproto            certain text proto file formats, see 'prototag' below
 -i ICTE=\<nn\>            ignore as many as \<nn\> connection table errors
 -i info                 collect extra text records in input (SDF and TDT only)
 -i skip=nn              skip over the first MM molecules in the input
 -i do=nn                only process NN molecules
 -i seek=offset          seek to byte offset OFFSET before starting reading
 -i stop=offset          stop reading once file is at byte offset OFFSET
 -i maxq=\<charge\>        set maximum plausible atomic partial charge
 -i minq=\<charge\>        set minimum plausible atomic partial charge
 -i mq=\<charge\>          set min (-charge) and max (+charge) plausible atomic partial charge
 -i mfc=\<charge\>         set min and max plausible formal charges
 -i d@3d                 Discern chirality from 3d coordinates
 -i ignore_bad_chiral    ignore obviously incorrect chirality specifications
 -i discard_chiral       discard all chiral information input
 -i mdlRincch            read incorrect chirality flags - explicit Hydrogen problem
 -i addmih               add implicit hydrogens to H deficient chiral centres
 -i dctb                 Discern Cis-Trans Bonds from geometry
 -i dwedge               Discern chirality from Wedge bonds only
 -i xctb                 Discard all cis-trans information on input
 -i ibctb                Ignore erroneous cis-trans bond input
 -i ignore_self_bonds    ignore bonds with the same atom at each end
 -i uccvno               immediately unconnect covalently bonded non organics on input
 -i Hiso                 allow implicit Hydrogens on aromatic isotopically labeled atoms [9c]
 -i rmhknown             remove implicit Hydrogens known property if possible
 -i RMHKNOWN             remove ALL implicit Hydrogens known attributes
 -i fixnd3v4             put formal charge on 3 connected, 4 bonded neutral Nitrogen
 -i mscale=\<scale\>       multiply all coordinates by \<scale\> upon reading
 -i mdlD                 allow D to be recognised as [2H]
 -i mdlT                 allow T to be recognised as [3H]
 -i mdl3chop             chop long elements in MDL files to first 2 chars
 -i mdl3=XX              change all long element symbols to XX
 -i ignore_bad_m         ignore unrecognised 'M' records in SDF files
 -i ignore_fatal_m       ignore otherwise fatal errors in M records
 -i mdlquiet             don't report unrecognised records to stdout
 -i SDFID:XX             name follows '\> \<XX\>' record in SDF file
 -i SDFNONAME            discard name in first record of SD file
 -i SDFNOPREPEND         do not prepend the SDF identifier to the identifier value
 -i SDFTAG2JSON          append all sdf tags as JSON form\n"
 -i firstsdftag          take first tag in an sd file as the name
 -i allsdfid             concatenate all sdf identifiers
 -i ISISEXTREG           try to discern the ISIS external registry number in SDF files
 -i mdlatomalias         replace elements by their alias symbols, A records in MDL files
 -i Galias               replace elements by their alias symbols, G records in MDL files
 -i gsubsdf=\<c\>          gsub all spaced values in sdf data records with \<c\>
 -i mdlRisonum           read isotopes as numbers rather than diffs from normal mass
 -i MDLIBT:num1=num2     translate input bond type 'num1' to 'num2'
 -i mdlustere            accumulate unclassified MDL chiral centres (type 3)
 -i mdlsep=\<..\>          separator between tags when reading mdl files 
 -i mdlnwce              do NOT write messages about chirality problems in sdf files
 -i sasge                MDL V30: convert single atom SGROUP labels to elements
 -i RDFID:XX             detect RDFILE info in tag \<XX\>
 -i RDFSTART:XX          records in RDFILES start with tag \<XX\>
 -i mol2fc               try to assign formal changes when reading MOL2 files
 -i mol2rfc              interpret mol2 charges as formal charges
 -i TDTSMI:XX            smiles is in tag XX rather than $SMI\<\>
 -i TDTID:XX             identifier is XX rather than PCN
 -i TDTAPPEND:XX         append dataitem XX to name
 -i TDTINCTAG            when appending TDT items, keep tags
 -i ignoretdtnosmi       ignore TDT's with no '$SMI' dataitem
 -i DOS                  remove '^M' characters from the end of each record
 -i delim=\<char\>         record delimiter (default '\n', use ^M for DOS)
 -i chemaxon             parse 'smiles |chemaxon| name' forms
 -i squoted              smiles records may have double quotes (typically "smiles |chemaxon|", name)
 -i csvsmicol=\<col\>      when reading csv files, smiles is in column number \<col\>
 -i csvsminam=\<col\>      when reading csv files, smiles is in column name \<col\>
 -i csvidcol=\<col\>       when reading csv files, identifier is in column number \<col\>
 -i csvidnam=\<col\>       when reading csv files, identifier is in column ame \<col\>
 -i csvallcol            when reading csv files, all extra columns become the name
 -i prototag=\<tag\>       when reading textptoto files, the smiles tag
 -i MOEFRAGS=\<fname\>     when reading Moe files, write named fragments to \<fname\>
```

## Input type
The most important aspect of reading molecules is to specify what kind of
file is being read. If your file ends with a recognised suffix, the tool
will assume the input type. But it is of course possible to do something
very strange such as
```
   -i smi file.sdf
```
which would read `file.sdf` line at a time, attempting smiles interpretation.
Don't do this!

If you are reading from stdin, by specifing `-` as the input file name, most 
tools will assume that the input is smiles. If it is not, then you will need to
specify what it is.

```
produce .smi output | tool -
produce .sdf output | tool -i sdf -
```
Multiple input types can be specified with one invocation
```
tool file.smi file.sdf file.gfp
```
is fine, but if you specify `-i` it will assume that specification
applies to **all** input files.

While there are a variety of input file types supported, some are
better supported than others. The types recognised are

* smi
* sdf 
* mdl  same as `sdf`
* rdf  reaction data format
* pdb
* mmod  macromodel file - ancient, might not work today.
* msi   old MSI/C2 file format. Obsolete.
* tdt   Thor DataTree - Daylight format.
* gfp   same as `tdt`
* mol2  Tripos mol2 format. Obsolete.
* psf
* crd
* mrk    Merck force field type
* chm    CharMM
* cif
* moe    not all features supported.
* mrv    not all features supported.
* inchi  only if compiled with conversion tools - not in public release.
* csv

The most heavily supported are `smi` and `sdf`. The others are really
best effort, but since they are not frequently used where LillyMol is
developed, they are probably not robust.

Specify the type of file being processed via `-i type` where `type` is one
of the items from the list above.

## Variations
Once the type of file being read has been established there are many variations
on behaviour that can be specified. There is a blurry line between a fix to the
connection table that is done as part of reading the molecule, and a fix
that should be done in the tool processing that molecule.

### -i ICTE
Ignore Connection Table Errors. This is the most frequently used `-i` qualifier.
As given above, **all** connection table errors in the input are ignored. Messages
will appear on stderr informing you of the problem. This way, you can process
a file containing 10 thousand molecules, all of which are invalid except
for one! This option should be used with caution.

If you want to limit the number of connection table errors tolerated, enter this
as `-i ICTE=nn` where a maximum of `nn` connection errors will be allowed,
and then a fatal error will happen.

Note that `fileconv` has a -B option what does the same thing as `-i ICTE`,
but it has the additional ability of being able to write the erroneous text
to a file `-B 999 -B log=bad.txt`. And of course `-i ICTE` also works for
`fileconv`.

### -i info
Some file formats contain extra records that are not part of the connection
table. MDL files are a comon example. Using the `-i info` option retains these
extra records as an array of strings in the Molecule. For most tools, 
they can be echo'd by using the `-o info` option, which writes that additional
text information to the output along with the connection table.

### -i skip=nn
Skip the first `nn` molecules in the file. These are skipped as text and
no Molecule interpretation is attempted. For a smiles file, individual records
are skipped, for a `sdf` file, data is skipped down the `$$$$` records.

### -i do=nn
Process only the first `nn` records in the file. Counting begins after any
`skip=` directive has been applied.

### -i seek=offset -i stop=offset
In an automated processign context, byte offsets into a file may have
been established, byte_offset_index, and we can direct that a seek operation happen
before reading begins. Note that if `offset` does not point to the 
start of a record, or more generally, to the beginning of a molecule,
unpredictable results should be expected.

Usually `seek=` is used on conjunction with `stop=` and when the
reader determines that it has reached the stop offset, it will terminate.

### Charges
Both formal and partial charges are subject to acceptable ranges, and
molecules containing values out of range are flagged as errors. The
default values can be adjusted.

### -i minq=\<charge\> -i maxq=\<charge\>
Separately set the min and max plausible partial atomic charge.

### -i mq=\<charge\>
Set both min and max partial chage to `-charge` and `+charge`.

### -i mfc=\<q\>
Set both min and max formal charge.

### -i dctb
If the molecule has a 2D or 3D geometry, use atom positions to discern
cis-trans bonding information. Note that LillyMol cannot canonicalize
cis-trans bonding. Perhaps that will be fixed one day.

### -i d@3d
Use the 3D coordinates of the input molecule to discern chiral centres.
Existing chiral centres specified are retained.

### -i ignore_bad_chiral
Ignore chirality errors that would otherwise be fatal.

### -i discard_chiral
Discard all chirality information as the molecule is read. Once the
instantiated Molecule is given to the tool, it contains no chirality.

### -i mdlRincch
Read incorrect chirality flags - explicit Hydrogen problem. This
seems not implemented. Do not use.

### -i addmih
Add implicit Hydrogen atoms to chiral centres that seem to be deficient. So
`C[C@](N)CC` would become `C[C@H](N)CC`. Although it should be noted that
the original representation is ambiguous. Over the years, a great many
smiles errors have been encountered.

### -i dwedge
Use wedge bonds to discern chirality. This generally works well if
there is just one wedge bond involved, but the algorithm does not
work if there are multiple wedge bonds to an atom. Basically I do
not know how to interpret that situation. Other software systems
yield varying results in this case. This is an open problem in
urgent need of information.

### -i ibctb
Ignore what might otherwise be classified as erroneous cis-trans bonding information.

### -i xctb
Remove all cis-trans bonding information in the input. Generally this is
a good idea. Even when present, such information is often not accurate.

### -i ignore_self_bonds
File formats like `mdl` can contain bonds with the same atom number
at each end. Ignore that bond, and attempt Molecule instantiation on
the remaining bonding records.

### -i uccvno
We have seen cases where acids might be represented as `O=C-[O-]-[Na+]`.
This option removes the bond between the Oxygen and the Sodium. More
generally, it removes any bond between a non-organic atom and an
organic atom. Note that this might have unintended consequences.

### -i Hiso
Add implicit Hydrogens to isotopic atoms that seem deficient. So
`C[9C]C` becomes `C[9CH2]C`'. This may avoid invalid valence rejections.

### -i rmhknown
In a smiles, remove the implicit Hydrogens known information if possible. So
`C[CH2]C` would be interpreted as `CCC`, but `C[CH]C` would not.

### -i RMHKNOWN
Unconditionally remove all implicit Hydrogens known information.

### -i fixnd3v4
Put a formal charge on 3 connected, 4 bonded neutral Nitrogen. This is a
common problem in connection tables, and can recover a molecule that would
otherwise be flagged as invalid.

### -i mscale=\<scale\>
Multiply all coordinates by \<scale\> upon reading. Useful for reading
data that might not be in Angstroms.

## MDL files
There are a great many options associated with reading MDL files.

### -i mdlD -i mdlT
Interpret `D` atoms as `[2H]`, Deuterium, and interpret `T` atoms as
`[3H]`, Tritium.

### -i mdl3chop
Chop long elements in MDL files to first 2 chars.

### -i mdl3=XX
Change all long element symbols to XX. Note that long symbols can also
be handled via `-E autocreate -E anylength`.

### -i ignore_bad_m
Ignore `M  ` records that are otherwise unrecognised. These are not
really bad, they are just unrecognised, the directive name is poorly
chosen.

### -i mdlquiet
By default, when unrecognised `M  ` records are encountered, that is
reported to stderr. That is seldom useful, and the `-i mdlquiet` option
suppresses that output.

### -i ignore_fatal_m
Ignore what might otherwise be fatal errors when processing `M  ` records.
Obviously, this is risky.

### -i SDFID:XX
The name of the molecule is in the `M  ` tag `XX`.

### -i SDFNOPREPEND
Normally when the name is extracted from a tag, the tag is prepended
to the name. So `-i SDFID:FOO` might produce a molecule with name
`FOO:bar`. The SDFNOPREPEND option results in the name `bar`.

### -i SDFTAG2JSON
Concatenate sdf tags into json form. So something that might be
entered as
```
>  <FOO>
bar
```
would be in the name as `{ "FOO": "bar" }`. Note that this option is largely
independent of other options, so if you want to collect all tags, you will
still need `-i allsdfid`, and you may also need `-i SDFNONAME`. Otherwise
the name field will be processed normally, leading to `smiles name { "FOO": "bar" }`
which might be exactly what you want.

### -I NAME2JSON
If needed, the molecule name can be encoded in json form as well as the
specified sdf tags. The resulting form will be
```
"name": "molecule name" "tag1": "something" ...
```

### -i SDFNONAME
Normally the first record in an MDL connection table is the name
of the molecule. If this is specified, discard that information instead.

### -i firstsdftag
The first tag in the MDL connection table is taken to be the name.

### -i allsdfid
The name is the concatenation of all tags in the MDL file.

### -i ISISEXTREG
Try to discern the ISIS external registry number in SDF files. Obsolete.

### -i IDM:TAG
This feature seems to have not been fully implemented. Ignore.

### -i RPSDFTAG=XX
replace the first sdf tag with \<XX\>. Seems to be no longer supported, do
not use.

### -i mdlatomalias
Replace elements by their alias symbols, A records in MDL files.

### -i Galias
Replace elements by their alias symbols, G records in MDL files.

### -i gsubsdf=\<c\>
Gsub all spaced values in sdf data records with \<c\>. Only applies to
data records coming after tags.

### -i mdlRisonum
Read isotopes as numbers rather than diffs from normal mass. Only applies
to V2000 files which have the isotope in the atom record.

### -i mdlsep=\<..\>
Separator between tags when reading mdl files. If a name is built up
from the concatentation of multiple tags, this is the separator between
those tags. Default is space.

### -i mdlnwce
Discerning chirality from wedge bonds and atom positions can often
be problematic. Tools will by default, issue warning messages about
these, but usually they are uninteresting. Suppress those messges.

### -i RDFID:XX
The name in an rdfile is in tag \<XX\>.

### -i RDFSTART:XX
Records in RDFILES start with tag \<XX\>.

### -i sasge
MDL V30: convert single atom SGROUP labels to elements.

## Mol2 files
These are obsolete files from Tripos.

### -i mol2fc
try to assign formal changes when reading MOL2 files.

### -i mol2rfc
Interpret mol2 charges as formal charges.

## TDT files
### -i TDTSMI:XX
smiles is in tag XX rather than $SMI\<\>

### -i TDTID:XX
identifier is XX rather than PCN

### -i TDTAPPEND:XX
Append dataitem XX to name

### -i TDTINCTAG
When appending TDT items, keep tags.

### -i ignoretdtnosmi
Ignore TDT's that do not contain a smiles. The next valid
TDT will be processed.

## CSV files.
Support for csv files is new and should be made more complete.

### -i csvsmicol=\<col\>
when reading csv files, smiles is in column number \<col\>.

### -i csvsminam=\<col\>
when reading csv files, smiles is in column name \<col\>.

### -i csvidcol=\<col\>
when reading csv files, identifier is in column number \<col\>.

### -i csvidnam=\<col\>
when reading csv files, identifier is in column ame \<col\>.

## Chemaxon files.

### -i chemaxon
Parse 'smiles |chemaxon| name' forms. Many, but not all Chemaxon extensions
are supported. This should be expanded.

### -i squoted
Smiles records may have double quotes (typically "smiles |chemaxon|", name).

## MOE files.
Support for Moe files is minimal.

### -i MOEFRAGS=\<fname\>
When reading Moe files, write named fragments to \<fname\>. This option is
problematic from a design perspective because it is not consistent with how
tools usually work. The problem is that there is no concept of separately
named fragments within the Molecule object. This is a kludge to get around
that. 

## Other options.

### -i DOS
remove '^M' characters from the end of each record. This is the default.

### -i delim=\<char\>
Record delimiter (default '\\n', use ^M for DOS) - historical only.

# Output
Within the LillyMol tools some tools will only generate smiles, while
others can generate a variaty of output forms. Generally we have tried
to reserve the `-S` option for specifying the name of an output file.

If the tool also supports a `-o` option, then it will have the ability to
write either a single output file, or multiple types. For example
```
fileconv ... -S newname -o smi -o sdf input.sdf
```
will read `input.sdf`, do whatever is specified on the command line
and write both `newname.smi` and `newname.sdf`.

Most tools with a -S option, will make an attempt to prevent overwriting
their input, although this is **NOT** uniformly enforced. So in some tools
it might be possible to do
```
tool ... -S foo -o smi foo.smi
```
which will first empty out `foo.smi`, and then start reading `foo.smi`,
which will happily report zero molecules in the input.

Tools which support a `-o` option should also support `-o help` which
yields
```
 -o smi3d       append coordinates after smiles
 -o flush       flush files after writing each molecule
 -o info        write any associated text info (if possible)
 -o MULT        create multiple output files, one per molecule
 -o MULT=nn     create multiple output files, \<nn\> molecules per file
 -o TOKEN=n     with MULT, use token N from the name as the file name
 -o ISIS        create as close as possible to an ISIS file
 -o AROM        when writing MDL files, write aromatic bonds as '4'
 -o V30         when writing MDL files, write V30 file type
 -o DOS         include ^M characters in output
 -o mdlWisonum  write isotopes as numbers rather than diffs from normal mass
 -o mdlWisoMnum write isotopes as numbers rather than diffs from normal mass
 -o mdlWincch   write stereo flags (1,2) without special consideration for
                explicit Hydrogens - may give incorrect values
 -o mdlMEND     always add the M  END record to MDL files
 -o RxSymbol    when writing R1 elements to MDL files, do NOT translate to R# form
 -o pdbso       write the atoms in PDB files in sequence (fragment) order
 -o pdbnws      label atoms in pdb files by number within sequence
 -o pdbec       label atoms in pdb files by element C1 C2 C3 N1 O1 O2 C4
 -o pdbsname    use stored atom names when writing PDB files
 -o mol2wfc     write the formal charge in the charge column of MOL2 files
 -o nochiral    exclude chirality info from smiles and mdl outputs
 -o nochiralflag don't write the chiral flag info to mdl files
 -o NOCT        exclude any CIS/TRANS information from smiles output
 -o smisep=\<c\>  separator between smiles and name: 'smisep=tab' for example
 -o \<type\>      specify one or more output types (smi,usmi,nausmi,rsmi,sdf,tdt,mol2,marvin)
```

Note that there is conceptual overlap between some `-o` options and `-i`
options. For example, `nochiral` will discard chirality just before the
molecule is written. This seldom makes sense, if you are not interested in
chirality, it should be discarded as the molecule is read, `-i nochiral`.

### -o smi3d
A LillyMol extension allows coordinates to be written in a smiles. This takes
the form `C{{{1.2,0.4,4.1}}`. Even though the option is named `smi3d` it works
for any coordinates. Some tools also support this via the `-K coords` directive.

### -o flush
Flush files after writing each molecule.

### -o info
Write any associated text info (if possible). This is usually paired with
the `-i info` flag, and only applies to certain file types, SDF and TDT.

### -o MULT
Write each structure to a separate file. Typically if a `-S` option is
specified, `-S foo` then `foo0.smi`, `foo1.smi` will be created, each 
containing a single molecule.

### -o MULT=nn
Multiple output files created, but this with with `nn` molecules in
each file.

### -o TOKEN=nn
When writing multiple files, normally the file name is the same
for all files created. If TOKEN=nn is specified, then the name
of each file created comes from a specific token in the name.
```
-o sdf -o smi -o MULT=10 -o TOKEN=2
```
would see output files created, each with 10 molecules. The name
of each file would be the second token in the name of each molecule
that started a group. And both `.sdf` and `.smi` files would be
created.

### -o ISIS
Create as close as possible to an ISIS file. The default `sdf` output
is an abbeviated form.

### -o AROM
When writing MDL files, write aromatic bonds as '4'.

### -o V30
When writing MDL files, write V30 file type always. By default that
only happens when the atom count exceeds 99.

### -o DOS
Include ^M characters in output.


### -o nochiralflag
Don't write the chiral flag info to mdl files. By default, if the molecule
contains chirality, the chiral flag is written.

### -o mdlWisonum
Write isotopes as numbers rather than diffs from normal mass

### -o mdlWisoMnum
write isotopes as numbers rather than diffs from normal mass

### -o mdlWincch
Write stereo flags (1,2) without special consideration for explicit Hydrogens -
may give incorrect values.

### -o mdlMEND
Always add the M  END record to MDL files.

### -o RxSymbol
When writing R1 elements to MDL files, do NOT translate to R# form.

### -o pdbso
Write the atoms in PDB files in sequence (fragment) order.

### -o pdbnws
Label atoms in pdb files by number within sequence.

### -o pdbec
Label atoms in pdb files by element C1 C2 C3 N1 O1 O2 C4.

### -o pdbsname
use stored atom names when writing PDB files.

### -o mol2wfc
write the formal charge in the charge column of MOL2 files.

### -o nochiral
Exclude chirality info from smiles and mdl outputs.

### -o nochiralflag
Do not write the chiral flag info to mdl files.

### -o NOCT
exclude any CIS/TRANS information from smiles output.

### -o smisep
By default, smiles are written as smiles followed by a space and then
the id. The space can be changed to any character with this directive.
```
-i smisep=vbar
```
will result in a vertical bar being the output separator.

### -o \<type\>
specify one or more output types (smi,usmi,nausmi,rsmi,sdf,tdt,mol2,marvin).
