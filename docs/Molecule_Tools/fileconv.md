# `fileconv`

## Background
`fileconv` reads, maybe changes, and maybe writes, streams of molecules.

`fileconv` is often the first part of a project, where things like
size constraints, undesirable atoms, and other kinds of filtering can
be applied. It is conceptually similar to a great many other similar
tools built on other libraries, obabel, molconvert, csfc are some that
spring immediately to mind.

`fileconv` clearly has too many options, and as a result, is unreasonably
complex. That said, every one of the options added to `fileconv` has been
added in order to meet a very specific, real-world need that has
come up.

Generally `fileconv` works well, but there are cases where the order
in which changes/checks are applied can be problematic. If that happens,
invocations can be pipelined.

## HOWTO
The **main** usage message is
```
Molecule_Tools/fileconv_opts.cc compiled  redacted redacted git hash 0608f40
usage: fileconv -i <input type> -o <output type> file1 file2...
The following options are recognised
  -f <directive> fragment selection, enter '-f help' for details
  -F <number>    exclude molecules having more than <number> fragments
  -O none        "organic" molecules only. '-O def' removes all bad elements
  -O <el>        "organic" molecules only, but 'el' is OK (repeat for each OK ele)
  -E <symbol>    create element <symbol>, use 'autocreate' for all
  -w <amw, -W <amw>  specify lower(-w) and upper (-W) amw limits
  -W LARGE       compute molecular weight based on largest fragment
  -c <number>    exclude molecules with atom count below <number>
  -C <number>    exclude molecules with atom count above <number>
  -C implicit    when computing atom count for -C, include implicit Hydrogens
  -r <n>         omit molecules with fewer than <n> rings
  -R <n>         omit molecules with more  than <n> rings. Use MRS:n for max rings in a system
                 both -r and -R take S:n for ring systems. Use 'spiro' to span spiro fusions
  -m <ringsize>  discard molecules with any ring larger than <ringsize>
  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed
  -Q <type>      compute partial charges - enter '-Q help' for info
  -I ...         isotope options, enter '-I help' for info'
  -s <qual>      chirality options, enter '-s help' for details
  -e             discard molecules with non periodic table elements
  -n ...         number assigner options, enter '-n help' for more info
  -o <type>      specify output file type(s), enter '-o help' for details
  -a             audit input only, no output
  -V             skip any molecule with abnormal valences
  -T <dx,dy,dz>  translate molecules (dx, dy, dz)
  -L <fname>     write rejected molecules to <fname>
  -B ...         handling for otherwise fatal problems, enter '-B help'
  -h <all>       make implicit hydrogens explicit
  -h <query>     make implicit hydrogens on atoms matching <query> explicit
  -S <string>    create output files with name stem <string>
  -p <string>    append various things to the name. Enter '-p help' for info
  -J <...>       fix obvious structure errors, enter '-J help' for info
  -Y <...>       miscellaneous options, enter '-Y help' for info
  -i <type>      specify input file type. Enter '-i help' for details
  -N <...>       Charge assigner specification, enter "-N help" for usage
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -H <..>        donor acceptor assignment, enter '-H help' for info
  -t E1=E2       standard element transformation options, enter '-t help'
  -A <qualifier> Aromaticity, enter "-A help" for options
  -K <...>       Enter "-K help" for SMILES options
  -v             verbose output
```

And several of those options expand to a large number of sub options.
The complexity of `fileconv` is overwhelming.

That said, most invocations call upon it to perform only a small subset
of operations possible. For example to get a quick idea of what might
be in a file of structures
```
fileconv -v -a rand.smi
2000 molecules had between 7 and 50 atoms. Average 28.0095
```
Is a common means of assessing a file of molecules.

## Flow
`fileconv` operates as a pipeline. Molecules are read from the 
input source. They are then processed, which may be a no-op, and
then, if requested, written to one or more output files.

`fileconv`, like all LillyMol tools, can read from stdin and write to
stdout, both of which are specified as `-`. So in order to read a
smiles from stdin, that might look like
```
something generating smiles | fileconv -v -a -
```
where the last `-` is interpreted as stdin.

This shows another important feature of `fileconv`. When it reads from
one or more files, it figures out what kind of file is being read
by examining the file name. Clearly `-` does not have a suffix,
so for this tool, smiles is assumed. If an sdf file is being
read from stdin, then the file type needs to be specified. For
example
```
something generating sdf | fileconv -v -a -i sdf -
```

`fileconv` output is specified with the `-S` option, and `fileconv`
can be part of a pipeline
```
generate smiles | fileconv <alter/filter smiles> -S - - | consume smiles
```
where both input and output (which defaults to smiles) are
from stdin and to stdout.

`fileconv` can also read and write gzip'd files
```
fileconv -v -a file.smi.gz
```
will work. The file type is discerned by stripping off the `.gz` and
looking at the suffix.

Like all LillyMol based tools, `fileconv` is designed to operate
silently. Adding one or more `-v` options will increase the verbosity
of logging messages.

By default, input files are taken from the command line, but if you
have a list of files in a file, you can provide that to fileconv via
```
fileconv ... F:/path/to/list/of/files
```
and it will read each file name to be processed from the `F:` file.

## Errors
All LillyMol tools, including `fileconv`, work on the assumption that
the input is correct and will die upon failing to read any input.
See the file `io.md` for more details on how `fileconv`, and all other
LillyMol tools behave.

It seems that the order in which the options are presented could be improved,
TODO:ianwatson. Investigate.

## Options

### -f ...
`fileconv` contains a number of directives for dealing with
fragment selection. When given `-f help` the following menu appears
```
  -f large       trim to largest fragment (can abbreviate to '-f l')
  -f alarge      determine largest fragment. Keep all fragments with that number of atoms
  -f lo          trim to fragment with most organic atoms
  -f lod         trim to fragment with most organic atoms and desirable features
  -f allo        keep all organic fragments
  -f RMDUP       remove duplicate fragments
  -f rmle=nn     remove fragments with NN or fewer atoms
  -f Q:qfile     keep largest  frag which matches query in <qfile>
  -f q:qfile     keep smallest frag which matches query in <qfile>
  -f Q:F:file    keep largest  frag which matches queries in <file>
  -f q:F:file    keep smallest frag which matches queries in <file>
  -f SMARTS:smt  keep largest  frag which matches smarts
  -f smarts:smt  keep smallest frag which matches smarts
  -f ALL:smt     keep all fragments that match smarts <smarts>
  -f rm:smt      remove all fragments that match <smarts>
  -f saltfile=<file> smiles file of known salts - always removed even if the largest fragment
  -f saltsmartsfile=<file> file containing smarts of known salts. All atoms in fragment must match (RDKit)
  -f parentfile=<file> file of known parent molecules - never removed as salts
  -f kmfok       compare known salts and parents by molecular formula only - not unique smiles
  -f kpallsalt   do not change a molecule if every fragment is a known salt
  -f rmxt=<n>    discard molecules with >1 fragment with more than n atoms
  -f rmxt        discard molecules with >1 fragment with more than 16 atoms
  -f sfs         sort fragments by size
  -f dmxt=<d>    discard molecules where largest fragments differ by <d> atoms or fewer
  -f manlf=<d>   discard molecules that have a non-largest fragment with more than <d> atoms
  -f klf=<d>     discard all but the <n> largest fragments
  -f RMF=<tag>   when processing TDT forms, write removed fragments to <tag>
  -f rmlarge     remove the largest fragment (arbitrary if two frags of same size
  -f rmlarge=<n> remove the largest <n> fragments (arbitrary if frags of the same size)
  -f rmsmall     remove the smallst fragment (arbitrary if two frags of same size
  -f rmsmall=<n> remove the smallst <n> fragments (arbitrary if frags of the same size)
  -f keepsmall   remove all but the smallest fragment (arbitrary if frags of the same size)
  -f keepsmall=<n> remove all but the smallest <n> fragments
  -f <number>    remove fragments so that all written molecules
                 have no more than <number> fragments
```
Such a menu of choices would normally be associated with a complex
stand-alone programme, but here it is a sub-menu of a single option
in fileconv. This is the story of fileconv. Over time, it has been
called upon to do more and more functions, and has become more and more
complex.

But as is often the case, it will seldom make sense to use more than 
one of these choices, although that is not a hard and fast rule.

### -f large
Trim to largest fragment (can abbreviate to '-f l'). Might not be the
best choice always.

### -f alarge
Determine largest fragment. Keep all fragments with that number of atoms.

### -f lo
Trim to fragment with most organic atoms. Prefer `lod` below.

### -f lod
Trim to fragment with most organic atoms and desirable features. This is
the current best choice for selecting what is likely to be the fragment
of interest. Like all heuristics it is not perfect.

Generally it will favour fragments with "organic" atoms, especially Nitrogen
atoms. So even though a fragment with just Carbon and Oxygen atoms might
have more atoms, it might instead select a fragment with a Nitrogen
atom as more likely to be the fragment of interest to a biological
study.

### -f allo
Keep all organic fragments. Any fragment that contains a non organic
atom will be discarded.

### -f RMDUP
Remove duplicate fragments. The unique smiles of each fragment will
be computed, and only unique forms retained.

### -f rmle=nn
Remove fragments with NN or fewer atoms. In complex mixtures might
be useful.

There are several query based fragment selection choices.

### -f Q:qfile     keep largest  frag which matches query in \<qfile\>
### -f q:qfile     keep smallest frag which matches query in \<qfile\>
### -f Q:F:file    keep largest  frag which matches queries in \<file\>
### -f q:F:file    keep smallest frag which matches queries in \<file\>
### -f SMARTS:smt  keep largest  frag which matches smarts
### -f smarts:smt  keep smallest frag which matches smarts
### -f ALL:smt     keep all fragments that match smarts \<smarts\>
### -f rm:smt      remove all fragments that match \<smarts\>

The syntax for specifying the query reflects what is available in
`tsubstructure`.

### -f saltfile=\<file\>
Smiles file of known salts - always removed even if the largest fragment.
This raises a fundamental question, how should counterions be handled?
Should there be a dictionary of known counterions, this option, or are
things better done via heuristics. Most use of fileconv relies on
heuristics, but there have been times when a saltfile has been needed.

###-f saltsmartsfile=\<file\>
Designed to read an RDKit salt file - those are smarts. They must match
all the atoms in the fragment. This makes no sense, but is included for
compatibility. Matching fragments are discarded.

### -f parentfile=\<file\>
File of known parent molecules - never removed as salts.

### -f kmfok
Compare known salts and parents by molecular formula only - not unique smiles.
This can be very efficient in some cases, but generally obscure use only.

### -f kpallsalt
Do not change a molecule if every fragment is a known salt.

### -f rmxt=\<n\> 
discard molecules with >1 fragment with more than `\<n\>` atoms. This is a 
very important option when dealing with molecules that are more likely
to be mixtures than single entities. I have found useful values for
this option to be numbers between 10 an 20. If a molecule contains
two fragments, each with more than (say) 15 heavy atoms, think of it
as a mixture and do not write.

### -f rmxt
Discard molecules with >1 fragment with more than 16 atoms. This just
provides a convenient default for `rmxt=\<n\>`.

### -f sfs
Sort fragments by size.  Does not do any fragment selection,
just re-orders the atoms on the molecule so that the fragments
appear sorted by size.

### -f dmxt=\<d\>
Discard molecules where largest fragments differ by \<d\> atoms or fewer.
This is another means of dealing with molecules that might reasonably
be thought of as mixtures.

### -f manlf=\<d\>
Discard molecules that have a non-largest fragment with more than \<d\> atoms.
Another way of dealing with molecules that might be mixtures.

### -f klf=\<d\>
Discard all but the \<n\> largest fragments.

### -f keepsmall=\<n\>
Remove all but the smallest \<n\> fragments.

### -f RMF=\<tag\> 
When processing TDT forms, write removed fragments to \<tag\>.

### -f rmlarge
Remove the largest fragment (arbitrary if two frags of same size).

### -f rmlarge=\<n\>
Remove the largest \<n\> fragments (arbitrary if frags of the same size).

### -f rmsmall
Remove the smallst fragment (arbitrary if two frags of same size).

### -f rmsmall=\<n\>
Remove the smallst \<n\> fragments (arbitrary if frags of the same size).

### -f keepsmall
Remove all but the smallest fragment (arbitrary if frags of the same size).

### -f \<number\>
Remove fragments so that all written molecules have no more than
\<number\> fragments. The first `\<number\> - 1` fragments are
removed, regardless of size.

Several different fragment selection criteria can be combined.

## -F 
The current `-F` option should be combined into the `-f` option
described above. It is seldom used and retained for compatability.
The usage message is
```
### -F \<number\>
Discard molecules with more than \<number\> fragments.

### -F mnsz=\<n\>
When counting fragments only count those with > \<n\> atoms.
```

## -O
Fundamental to LillyMol is the idea of Organic elements. By default this set
includes H, C, N, O, F, P, S, Cl, B, I. Using `-O def` will remove
all molecules that, after fragment processing, still contain non organic
elements. That way molecules that contain a Sodium salt, can be
selected for as organic. Without fragment selection, that molecule would
be rejected for containing non organic atoms.

Some elements are sometimes thought of as Organic, Selenium for example,
Boron is another. If you want to have a given element considered organic
during an invocation, do this with the `-O Se` option.

Note that this can also be accomplished via the `-E` option in almost
all tools. `-E O:Se` works for all tools. The `-O` option was added
to `fileconv`, and once the utility was recognised, adding the
functionality to `-E` made it more available. See discussion of
elements.

Boron is unusual in that some forms of Boron may be considered OK
in an "organic" subset of molecules. The `-O okOBO` option combination
specifies that Boron atoms are Ok if they are in the form 'O-B-O. All
other Boron atom types are classified as non organic.

Generally fragment stripping directives are processed before organic
element discovery, so '-f lod -O none' will pass '[U].CC' as 'organic'.

## -E ...
Allow specification of elements. See separate doc.

## -w \<amw\> -W \<amw\>
Specify lower (-w) and upper (-W) molecular weight limits. Molecules
falling outside these ranges are discarded. Note that the calculation
is performed after fragment selection.

## -W LARGE
Compute molecular weight based on largest fragment. Make sure there
is no fragment selection active.

### -c \<number\>    exclude molecules with atom count below \<number\>.
### -C \<number\>    exclude molecules with atom count above \<number\>.
With both options, only explicit atoms are considered, so beware
if your molecules contain explicit hydrogens. The filtering is
done after any operations that might remove explicit hydrogens are
performed.

### -C implicit
When computing atom count for -C, include implicit Hydrogens.

### -r <n> omit molecules with fewer than <n> rings
### -R <n>  omit molecules with more  than <n> rings.
Use MRS:n for max rings in a system.
Both the `-r` and -R take S:n for ring systems.
Use 'spiro' to span spiro fusions.

```
filtering -r 2 -R 4
```
will discard molecules that have fewer than 2 or more that 4 rings.
```
fileconv -r S:2 -R S:4
```
will discard molecules that have fewer than 2 or more than 4 ring
systems.
```
fileconv -r S:2 -R S:4 -r spiro
```
extends the definition of a ring system to spiro fused rings.

Lower limits on the number of aliphatic and aromatic rings can
be imposed.
```
fileconv -r aliph:1 -r arom:2 ...
```
will only write a molecule if it has at least 1 aliphatic ring and
at least two aromatic rings. The same syntax applies for the `-R`
option, except that then the values apply to maximum numbers of
aliphatic and aromatic rings.

### -m \<ringsize\>
Discard molecules with any ring larger than \<ringsize\>. Note that
like many filtering options in `fileconv` the same thing could
be done in `tsubstructure`.

### -X \<symbol\>
Extract/remove all atoms of type \<symbol\>. No bonds changed.
Most frequently this is used for removing explicit Hydrogens. Note
however that explicit Hydrogens are also removed in chemical
standardisation.

### -Q \<type\>
Compute partial charges - enter '-Q help' for info. LillyMol
has several partial charge calculations built in. You will need
to be writing to a file format that supports partial charge
information. This is infrequently used.

### -I Isotopes.
Over the years, isotopes have proven to be very convenient atom
markers. In LillyMol, an isotope is simply a positive integer
with no relation to any concept of protons and neutrons. As such
a great many operations can be performed on isotopes. Entering
`-I help` yields the following options.
```
 -I 0             discard molecules containing any isotopic atoms
 -I change        change any isotopic atoms to normal form
 -I change=<n>    change any isotope <n> atoms to normal form
 -I change=<i,j>  change any isotope <i> atoms to isotope <j>
 -I CHANGE        change to normal form. Free implicit H. Should be default
 -I alliso=<i>    change any isotopic atoms to isotope <i>
 -I smarts:smarts,<i>    change any atoms matching <smarts> to isotope <i>, e.g. '[2C],0' or '[#16],4'
 -I firstatom:smarts:<i> in a multi-atom smarts, set the isotope only on the first matched atom
 -I remove        remove any isotopically labelled atoms
 -I add:\<i\>-\<smiles\> attach a fragment to all occurrences of isotope \<i\>

```

### -I 0
Discard molecules containing any isotopic atoms. Note that an atom
like [12CH4] is considered isotopic, even though 12 is the normal
isotope for Carbon. If the atom has the isotope attribute set, it is
isotopic. `-I none` is also recognised.

### -I change
Change any isotopic atoms to normal form. Note that this does
not necessarily clear the implicit Hydrogen known flag. The
mroe recent `CHANGE` directove below should be preferred. The
only reason it is not the default is for historical reasons.

### -I CHANGE
Change to normal form. Free implicit H. Should be default.

### -I change=\<n\>
Change any isotope \<n\> atoms to normal form. This allows treating
different subsets of atoms in different ways. Note that `convert` is
also recognised as a synonm for `change`.

### -I change=\<i,j\>
Change any isotope \<i\> atoms to isotope \<\j>. A more natural
syntax would have been to use `=` instead of `,`. TODO:ianwatson
enable that.

### -I alliso=\<i\>
Change any isotopic atoms to isotope \<i\>.

### -I smarts:smarts,\<i\>
Change any atoms matching \<smarts\> to isotope \<i\>, 
e.g. `[2C],0` or `[#16],4`.

### -I remove
Remove any isotopically labelled atoms. May result in
multiple fragments in the resulting molecule.

### -I add:\<i\>-\<smiles\>
A common operation is to need to append a fragment to an 
isotopically labelled molecule - usually an isotopically labelled
reagent. By convention, the first atom in the smiles is the
attachment point. So if you had a set of acids where the oxygen
had been removed and isotope 12 left on each of the remaining
carbon atoms, reconstruct the acid with `-I add:12-O`.


## -s Steroechemistry
Note that this is another example where certain functionality was
initially implemented in `fileconv`, but some of that functionality
later migrated to the `-i` option, and is available to all tools. See
the `io` file.

Entering `fileconv -s help` yields
```
  -s good        ignore erroneous chiral input
  -s find        identify all stereo centres in the molecules (random chirality)
  -s rfind       identify stereo centres for ring atoms only
  -s invert      invert all stereo centres
  -s 1           include chiral info in output (default)
  -s 0           exclude chiral info from output
  -s remove      remove all chirality upon read
  -s rmchiral=... remove chiral centres on atoms that match query(s)
  -s 'rmchiral=SMARTS:[ND2H]' will remove 2 connected Nitrogen chiral centres
  -s 'rmchiral=SMARTS:[AR3]' will remove wedge bonds in adamantyl systems
  -s rmbad       remove chiral centres that really aren't chiral
  -s xctb        remove all cis-trans bonding information
  -s rmbctb      remove invalid cis-trans bonding information
  -s nonocm      no non organic chirality messages
  -s rmnoc       remove chirality from non organic atoms
  -s maxc=<n>    discard molecules with more than <n> chiral centres
```

### -s good
Ignore erroneous chiral input.

### -s find
Identify all stereo centres in the molecules (random chirality). There
is a library function, `is_actually_chiral` that checks whether or not
an atom is truly asymmetric. If it is asymmetric, and it is not marked,
then a random chirality will be assigned. This is more for testing 
purposed.  See `id_chirality` for a tool that enumerates unlabelled
chiral atoms.

### -s rfind
Identify stereo centres for ring atoms only.

### -s invert
Invert all stereo centres.

### -s 1
Include chiral info in output (default).

### -s 0
Exclude chiral info from output.

### -s remove
Remove all chirality upon read. If chirality will be removed on
output, you may as well use `-s 0` to remove it during input.

### -s rmchiral=...
Remove chiral centres on atoms that match query(s).

```
-s 'rmchiral=SMARTS:[ND2H]'
```
Will remove 2 connected Nitrogen chiral centres

```
-s 'rmchiral=SMARTS:[AR3]'
```
Will remove wedge bonds in adamantyl systems.

### -s rmbad
Remove chiral centres that really aren't chiral. This also uses
`is_actually_chiral` to determine if labelled atoms are actually
asymmetric.

### -s xctb
Remove all cis-trans bonding information.

### -s rmbctb
Remove invalid cis-trans bonding information.

### -s nonocm
No non organic chirality messages. By default, when chirality is
found associated wth non organic atoms, a warning message is issued.
This suppresses that message.

### -s rmnoc
Remove chirality from non organic atoms.

### -s maxc=\<n\>
Discard molecules with more than \<n\> chiral centres.

## -e
Discard molecules with non periodic table elements. Note that if
`-i mdlD` and/or `-i mdlT` are in use, then those atoms are interpreted
as Hydrogen isotopes, and so are in the periodic table.

## -n ...
Number assigner. If you have molecules that do not have an id, this is
a means of generating an id for them.

The usage message is
```
  -n <number>    assign sequential numbers R(%d) starting with <number>
  -n <string>    use <string> rather than 'R'
  -n 'replace'   replace the entire name with R(nnn)
  -n xrf=fname   create a cross reference file
  -n noparen     omit parentheses when forming the names
  -n digits=<n>  left pad numbers with leading 0's to yield <n> digits
  -n nochange    only apply if the existing name is empty
  -n nsep='c'    separator between number and existing name (default ' ')
  -n seq         assign unadorned sequential numbers - no prefix or parentheses
```
This is seldom used, and is not described in detail here.

The most simple use case would be for a set of molecules that did not contain
identifiers, or had duplicate identifiers, and you need to apply a unique
identifier to each molecule.
```
fileconv -n seq ....
```
will apply sequential numbers to the molecules, inserting a sequential
number at the beginning of each record.
```
fileconv -n seq -n replace ...
```
will discard the name entirely and only the sequential numbers will appear.

Beware 'naked' numbers which may get confused with other numbers. The default
for the number assigner is to create ids that look like 'R(nnn)'.

## -o
The `-o` option governs what kind of output is produced. This is covered
in the `io` document.

## -a
Audit input only, no output. This is useful for getting a quick idea
of the atom counts in a file.
```
fileconv -v -a file.smi
```
tends to be a common operation. This will also let you know if the
smiles are valid using default options.

## -V
Skip any molecule with abnormal valences. Note that this is subject to
the data in the in-built periodic table. There is not much flexibility
in altering that. Note too that most elements do not have valence
information, and so cannot fail valence tests. So
```
[U](C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)(C)C
```
will pass the -V option. Note however the maximum allowable
formal charge concepts defined int he `io` document.

## -T \<dx,dy,dz\>
Translate molecules (dx, dy, dz). Only works if the starting
molecule has coordinates.

## -L \<fname\> 
Write rejected molecules to \<fname\>. For historical reasons
this will create two files, one containg just ids and the
other a structure file of the of rejected structures. The
file type will be as specified by the `-o` option if present,
smiles otherwise.

## -B ...
Handling for otherwise fatal problems, enter '-B help'.
After the usefullness of being able to skip otherwise
fatal errors was realized, the `-i ICTE` option was added.
The `-B` option remains in fileconv for historical reasons, 
although it does have some extra functionality.

The usage message is
```
 -B <nn>        ignore as many as <nn> otherwise fatal input errors
 -B log=<fname> echo (unchanged) rejected input records to <fname>
```

The first can be replaced by `-i ICTE=\<nn\>`, but the second
is unique to `fileconv`. As bad connection tables are encountered
it will seek back to the start of that bad connection table, and
echo it, as text to the file designated. When dealing with
large files this can be very helpful in tracking down problems.

## Implicit Hydrogns to Explicit Hydrogens.
The `-h` option controls this

### -h \<all\>
Make implicit hydrogens explicit.
### -h \<query\>
Make implicit hydrogens on atoms matching \<query\> explicit. THis
follows standard query specifications, `-h SMARTS':[CR0]'` would
only make implicit hydrogens explicit if they are on non ring carbon
atoms. By default, it is an old style substructure query file.
### -h last=\<z\>
Move all atoms with atomic number `z` to the end of the connection
table. This is often useful when writing .sdf files where putting
the hydrogen atoms last in the atom list is frequently seen.

## -S \<fname\>
Specify the file name stem for any output. For each output type
specified with the `-o` option, a corresponding `\<fname\>.type` file
will be created.

## Properties
It has proven convenient to be able to do certain molecular 
calculations with fileconv. This is not a great idea, since `iwdescr` is
the proper tool for that. But the convenience of `fileconv` has been
very compelling. 

Entering `fileconv -p help` yields
```
  -p AMW         append molecular weight to the name (isotopic atoms fail molecule)
  -p AMWI        append molecular weight to the name (isotopic atoms handled)
  -p MF          append molecular formula to the name
  -p ISISMF      append ISIS-like molecular formula to the name
  -p NATOMS      append number of atoms to the name
  -p NRINGS      append number of rings to the name
  -p AROMR       append number of aromatic rings to the name
  -p HTROAC      append number of heteroatoms to the name
  -p EXACT       append exact mass to the name
  -p NFC         append net formal charge to the name
  -p CLND        append Chemiluminescent Nitrogen Detection
  -p PREPEND     do a prepend rather than append
  -p LARGE       computed properties derived from the largest fragment
```
These are not described here, but an output might look like
```
C(=O)NCCCC CHEMBL45466 NATOMS = 7 NRINGS = 0 AMW = 101.147

```
Again, `iwdescr` is a more comprehensive and better formatted tool for
getting molecular properties.

## Fixing obvious structure errors, -J
The following options are available under the `-J` option.
```
  -J nitro     changes -N(=O)-O to a Nitro
  -J pyridine  puts a + charge on substited Pyridines
  -J pquatn    puts a + charge neutral Quaternary Nitrogens
  -J SD3+      puts a + charge on 3 connected, saturated Sulphur
  -J XD2+      puts a + charge on 2 connected, saturated Halogens (Cl, Br, I)
  -J app=xxx   append xxx to each changed structure
  -J rmhsqb    remove implicit hydrogens known flag from otherwise OK atoms
  -J rmbcharge try removing charges on atoms with bad valences
  -J Xnc       remove formal charges from halogen atoms
  -J Cnc       remove formal charges from carbon atoms
  -J all       turn on all known changes
```
These transformations represent structural errors that have occurred
multiple times over the years, and it has been worth the effort to fix
those molecules, rather than being forced to discard them. There
have of course been a great many other structural errors encountered,
but these are the only ones for which we have chosen to implement
a systematic fix.

## Misc
As `fileconv` was called upon to do more and more things, we quickly
ran out of letters in the alphabet to use for options, so a great deal
of functionality is available behind the `-Y option. Here is the
usage message
```
 -Y nbvm          No Bad Valence Messages, or strange electron counts
 -Y okbvi         during valence check, ok to have bad valence on isotopes
 -Y appchg=xxxx   append 'xxxx' to the name of molecules that are changed
 -Y pblen         print all bond lengths in the molecules
 -Y pbang         print all bond angles in the molecules
 -Y ptor          print all torsion angles in the molecules
 -Y pmaxd         print the max interatomic distance in each molecule
 -Y dbg           debug print each molecule
 -Y namerx=<rx>   discard molecules unless the molecule name matches <rx>
 -Y grepvname=<rx> discard molecules if the molecule name matches <rx>
 -Y ftn           keep only the first token in molecule names
 -Y nsubw=c       translate all whitespace in molecule names to 'c'
 -Y chname=rx     change name to what is matched by <rx> (optional match)
 -Y CHNAME=rx     change name to what is matched by <rx> (MUST match)
 -Y tfirst=char   truncate name at first <char>
 -Y tlast=char    truncate name at last <char>
 -Y NPREPEND=s    prepend <s> to each name
 -Y ntoken=n      the output name is word 'n' in the input name
 -Y maxplen=<n>   discard molecules with max path length > <n>
 -Y aclf          atom counts are for the largest fragment only
 -Y nhsqb         explicit hydrogen atoms in smiles written without square brackets
 -Y rmsqb         remove unnecessary square bracketed atoms  - Hcount is OK as specified
 -Y FHrmsqb       in atoms like [C] free the H count to what is computed. Square brackets removed
 -Y Xihaltvalence discard molecules like =S- where an implicit hydrogen satisfies an alternate valence
 -Y rmamap        remove atom map numbers
 -Y fixarom=smarts call find_kekule_form on the ring (system) matched by the first atom in <smarts>
 -Y rmatom=smarts remove all atoms that match <smarts>
 -Y zpad=width    left pad the name with 0's to <width>. Negative numbers remove leading chars from name
```
### -Y nbvm
No Bad Valence or strange electron counts warning messages.

### -Y okbvi
During valence check, ok to have bad valence on isotopes.

### -Y appchg=xxxx
Append 'xxxx' to the name of molecules that are changed. Molecules
passing through `fileconv` can be changed either by `fileconv` itself
or in chemical standardisation.

### -Y pblen
Print all bond lengths in the molecules. A debugging aid.

### -Y pbang
Print all bond angles in the molecules. A debugging aid.

### -Y ptor
Print all torsion angles in the molecules. A debugging aid.

### -Y pmaxd
Print the max interatomic distance in each molecule. A debugging aid.

### -Y dbg
Debug print each molecule. A debugging aid. Voluminous output may result.

### -Y namerx=\<rx\>
Discard molecules unless the molecule name matches \<rx\>. If you
have smiles, just use grep. There should be functionality to query
the associated information sometimes present in .sdf files, but that
does not exist today.

### -Y grepvname=\<rx\>
Discard molecules if the molecule name matches \<rx\>. If you have
smiles, just use grep.

### -Y ftn
Keep only the first token in molecule names.

### -Y nsubw=c
Translate all whitespace in molecule names to 'c'.

### -Y chname=rx
change name to what is matched by \<rx\> (optional match). There
must be a single capture group in the regular expression.

For example, given 
```
cat /tmp/foobar.smi
C fooQQQbar
fileconv -Y 'chname=^foo(..+)bar$' -S - /tmp/foobar.smi
```
yields `C QQQ`. If you have smiles, `sed` will likely be easier.

generate

### -Y CHNAME=rx
change name to what is matched by \<rx\> (MUST match).

### -Y tfirst=char
Truncate name at first \<char\>. If using smiles, `sed` is likely
easier.

### -Y tlast=char
Truncate name at last \<char\>.

### -Y NPREPEND=s
Prepend <s> to each name. If using smiles, `sed` is likely easier.

### -Y ntoken=n
The output name is word 'n' in the input name. If using smiles, `cut`
may be easier to use.

### -Y maxplen=\<n\>
Discard molecules with max path length > \<n\>. Note that filtering
by any computed property is available in `iwdescr`.

### -Y aclf
Atom counts are for the largest fragment only.

### -Y nhsqb
Explicit hydrogen atoms in smiles written without square brackets.
We have observed different requirements in different smiles implementations.

### -Y rmsqb
Remove unnecessary square bracketed atoms  - Hcount is OK as specified.
Given `[CH4]` that will be written as `C` since the `H4` satisfies the
valence and the atom is unambiguous without the brackets.

### -Y FHrmsqb
In atoms like [C] free the H count to what is computed. Square brackets removed.
That atom will be written an `C`. The assumption is that the Hydrogen count
in the starting smiles was omitted.

### -Y Xihaltvalence
Discard molecules like =S- where an implicit hydrogen satisfies an
alternate valence. Some elements have alternate valence, with Sulphur being
the main culrit. In the case above, there is al alternate valence of 4, and
so an implicit Hydrogen would be placed on the `S` atom. In this case,
consider this a valence error instead - since it is satisfying an
alternate, rather than normal valence.

### -Y rmamap
Remove atom map numbers.

### -Y fixarom=smarts
Call find_kekule_form on the ring (system) matched by the
first atom in \<smarts\>. The matched atoms are found. Then the ring system
containing those atoms are identified. All bonds in those ring systems are
set to single bonds, and an aromaticity determination is attempted. This
helps where an erroneous bonding specification has prevented a ring from
being perceived as aromatic.

### -Y rmatom=smarts
Remove all atoms that match \<smarts\>. This may result in multiple
fragments being formed. Use `trxn` for more control.

### -Y zpad=width
Left pad the name with 0's to \<width\>. So given `C 123` and
`fileconv -Y zpad=4` the result will be `C 0123`.

Negative numbers remove leading chars from name, so given `C X123` and
`fileconv -Y zpad=-1` the result will be `C 123`. If using smiles,
`sed` may be easier.

### -t E1=E2
Element transformations. For example, `-t I=Cl -t Br=Cl` transforms
all heavy halogens to Chlorine.

## -i \<type\>
Specify the input type. See discussion in `io`.

## -N ...
Formal charge assigner. Assigns formal charges based on rules
developed by Fred Bruns.

## -H ...
Hydrogen bond donor/acceptor assignments based on rules from
Fred Bruns.

## -g Chemical Standardisation
See separate document on chemical standardisation.

## -A Aromaticity.
See separate document on aromaticity.

## -K Smiles options.
Several tools support a `-K` option which enables control over
various aspects of smiles parsing and formation. The usage
message is.
```
  -K nru         do not reuse ring closure numbers when outputting smiles
  -K ctb         include cis-trans bonds in the smiles
  -K <number>    use <number> as a seed for random smiles
  -K nochiral    exclude chiral info from all smiles produced
  -K coords      include coordinates in the smiles - '{{xx,yy,zz}}'
  -K Hiso        provide implicit Hydrogens to isotopic atoms needing H atoms
  -K nhis        exclude implicit Hydrogen information from smiles
  -K iusmi       include isotopic information when determining canonical order
  -K niusmi      exclude isotopic information when determining canonical order
  -K dbusmi      include directional bonding information when determining canonical order
  -K ndbusmi     exclude directional bonding information when determining canonical order
  -K rnoff=nn    add NN to all ring numbers
  -K nonH        do NOT write implicit Hydrogen info on pyrrole-like N atoms
  -K nihc        do NOT consider implicit Hydrogens during canonicalisation
  -K esssr       perceive the extended SSSR ring set
  -K noD         do NOT include the D operator in smarts generated
  -K iso01       when computing unique smiles, isotopes are zero or non zero
  -K rcsbd       include directionality in ring closure single bonds
  -K nd4h        four connected neutral Nitrogen atoms have a Hydrogen
  -K fcnum       write formal charges as +[number] rather than ++
```

A repeated theme here is that there are often several ways of doing
the same thing. These options have their effect either when molecules
are built from smiles, or when a molecule is producing a smiles.

### -K nru
Do not reuse ring closure numbers when outputting smiles.

### -K ctb
Include cis-trans bonds in the smiles. Note that cis trans bonding
information can also be controlled with the `-i` option.

### -K \<number\>
Use <\number>\ as a seed for random smiles. If random smiles are
written, this controls the smiles order.

### -K nochiral
Exclude chiral info from all smiles produced. Note the interaction
with the `-i` option which can exclude chirality from the input.

### -K coords
Include coordinates in the smiles - '{{xx,yy,zz}}'.

### -K Hiso
Provide implicit Hydrogens to isotopic atoms needing H atoms. Note
the interaction with some of the `-Y` options to `fileconv`, but the
`-K` option will be available on other tools that make a `-K` option
available.

### -K nhis
Exclude implicit Hydrogen information from smiles. All atoms that can be
written without square brackets will not have them.

### -K iusmi
Include isotopic information when determining canonical order. This
is on by default.

### -K niusmi
Exclude isotopic information when determining canonical order. Normally
isotopic information is used to form the canonical atom ordering. Use of
this may lead to non unique, unique smiles.

### -K dbusmi
Include directional bonding information when determining canonical order.
Directional bonding in LillyMol does not work very well, and generally
this should not be used.

### -K ndbusmi
Exclude directional bonding information when determining canonical order.
This is recommended.

### -K rnoff=nn
Add NN to all ring numbers. Mostly for testing and debugging.

### -K nonH
Do NOT write implicit Hydrogen info on pyrrole-like N atoms. Different
smiles implementation may have different requirements, or we may
deliberately wish to introduce ambiguity.

### -K nihc
Do NOT consider implicit Hydrogens during canonicalisation. Using this
may result in unique smiles failure.

### -K esssr
Perceive the extended SSSR ring set. Testing and debugging only.
This will result in rings that would otherwise be discarded as
redundant being found. It is not really specific to smiles formation.

### -K noD
Do NOT include the D operator in smarts generated. A LillyMol Molecule
can produce a smarts. Such a smarts is made useful by **not** specifying
the D qualifier in the atomic smarts.

### -K iso01
When computing unique smiles, isotopes are zero or non zero. Mostly
for testing and debugging.

### -K rcsbd
Include directionality in ring closure single bonds.

### -K nd4h
Four connected neutral Nitrogen atoms have a Hydrogen.

### -K fcnum
write formal charges as +[number] rather than ++.
