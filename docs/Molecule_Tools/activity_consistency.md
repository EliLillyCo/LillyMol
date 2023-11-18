# Activity Consistency

# Background

Many QSAR datasets contain duplicate data, either
exact duplicates, or when counterions or chirality is discarded. The
task of `activity_consistency` is to examine the nature of those
inconsistencies, and optionally create a new activity file that
contains a composite view of the data.

# HOWTO
The following options are recognised.
```
Groups identical molecules and compares activities
  -X <fname>    experimental data
  -e <col>      activity is column <col> in the -X file
  -e <n>        activity is token <n> in the molecule name (no -X)
  -a            reduce to graph form
  -c            ignore chirality
  -x            ignore cis-trans bonds
  -m            only display common structure record of groups
  -t <tol>      activity difference for items to be merged
  -r            interpret the -t option as a ratio between activities
  -O <fname>    write groups outside tolerance to <fname>
  -N <fname>    write non duplicates to <fname> - cannot use with -t,-O
  -M <stem>     write new smiles and activity files with max activity 
  -M max        write max activity to the -M file
  -M ave        write average activity to the -M file
  -M median     write median activity to the -M file
  -M rand       write a random example to the -M file
  -M range      write a random value from within the range to the -M file
  -M rmdup      discard all duplicate structures
  -M rminc      discard structures with inconsistent class assignments
  -T ...        standard element transformation options
  -V ...        what to do with multi-valued activities (-V help for info)
  -K <ele>      remove elements of type <ele>
  -W ...        specification of allowed elements
  -U <fname>    tabular output format
  -b <natoms>   lower atom count cutoff
  -B <natoms>   lower atom count cutoff
  -z            remove leading 0's from identifiers
  -l            reduce to largest fragment
  -i <type>     input specification
  -g ...        chemical standardisation options
  -E ...        standard element specifications
  -A ...        standard aromaticity specifications
  -v            verbose output
```
The activity files are all assumed to be space separated - that should
be made more flexible.

The experimental data can be specified either via the `-X` option, 
or if it is a particular column in the smiles file, via the `-e` option
with no `-X` option.

## TLDR
A common use of this tool might look like
```
activity_consistency -E train.activity -e 2 -l -c -x -g all -T I=Cl -T Br=Cl
        -M max -M merged -v train.smi
```
New files `merged.activity` and `merged.smi` are created, containing
merged information about structural duplicates.

## Options

### -a
Reduce to graph form. All bonds single, all atoms carbon. More likely
useful for exploration only.

### -c
Ignore chirality. Since chirality information is often missing or
inaccurate, this is usually recommended.

### -x
Ignore cis-trans bonds. Especially since LillyMol does not properly
canonicalize cis-trans bonds. This should always be used.

### -m
Only display common structure record of groups. Writes a summary of
the id's grouped together, rather than listing them all separately.

### -t \<tol\>
Activity difference for items to be merged. The `-O` option contains
all structure groups that are designated as 'inconsistent'. In order
to do that with a continuous response, we need to specify how
different do values need to be in order to be considered 'different'.

### -r
Interpret the `-t` option as a ratio between activities.

### -O \<fname\>
Write groups outside tolerance to \<fname\>. See the `-t` and `-r`
options for specifying the tolerance.

### -N \<fname\>
Write non duplicates to \<fname\> - cannot use with `-t` , `-O`.

### -M \<stem\>
Write new smiles and activity files with max activity. This is the
most common usage for `activity_consistency`. The following qualifiers
specify what gets written to the `-M` file.

### -M max
Write max activity to the -M file. This usually seems to be the best
idea. For example, where chirality might be unknown, Chemists will
figure out which enentiomer is the active one and will use it.

### -M ave
Write average activity to the -M file.

### -M median
Write median activity to the -M file.

### -M rand
Write a random example to the -M file.

### -M range
Write a random value from within the range to the -M file.

### -M rmdup
Discard all duplicate structures.

### -M rminc
Discard structures with inconsistent class assignments.

### -T ...
Standard element transformation options. Most commonly this is to
group the heavy halogens `-T 'I=Cl' -T 'Br=Cl'`. Differentiating
among these seldom leads to a better model.

### -V ...
What to do with multi-valued activities (-V help for info). The
activity file may contain duplicate data - same identifier, possibly
different values.
```
12345 3.14
12345 -3.14
```
or different class assignments
```
12345 foo
12345 bar
```
The following qualifiers are allowed

#### -V first
Take the first instance of data.

#### -V ave
Take the average activity where multiple instances exist.

#### -V max
Take the maximum activity where multiple instances exist.

#### -V WRITE=\<fname\>
Write multiple instances to \<fname\>.

#### -V append
When writing a file containing duplicate identifiers, concatenate
the identifiers.

### -K \<ele\>
Remove elements of type \<ele\>. This is an obscure option. It might be
useful for removing explicit Hydrogens, but those can also be removed
using chemical standardisation, `-g all`.

### -W ...
Specification of allowed elements. By default, only organic elements
are permitted in the output. Other elements can be temporarily considered
as organic, `-W Si`.

### -U \<fname\> 
For all of the equivalent structure groups, writes a
tab separated record showing the folling information about that 
set of molecules. The `ID` will be that of the first molecule
encountered within that group.
```
Minval Maxval Range Average ID Activity
```

### -b \<natoms\>
Lower atom count cutoff. Discard molecules having too few atoms. Many
datasets will contain molecules that are either 'very small' or 'very large'.
These outliers may have negative influences on model performance, and we
usually find it useful to suppress them.

### -B \<natoms\>
Upper atom count cutoff.

### -z
Remove leading 0's from identifiers. This may be useful in an environment
where sometimes identifiers are present as '000000123456' and sometimes as '123456'.

### -l
reduce to largest fragment.
