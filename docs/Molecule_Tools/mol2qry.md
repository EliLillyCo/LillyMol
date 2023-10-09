# mol2qry
# Purpose
`mol2qry` converts molecules to substucture query forms.

A frequent task is the need to take structural information
observed somewhere, and use that to interrogate a larger collection
of molecules. Sometimes this can be done via a similarity search,
possibly using Tversky parameters, but other times an exact
match is required. Especially if there are restrictions on how
the smaller molecule is embedded in the larger molecule.

Another very common workflow is to have a set to molecules that
can only be substituted at a given atom, something that could not
be done with a similarity search.

# HOWTO
The usage message is.
```
Molecule_Tools/mol2qry.cc compiled redacted redacted git hash ec071c7
Converts a molecule to a query file
  -m             all ncon and nbonds values are written as minima
  -r             all ring bonds become type ANY
  -j             all atoms conserve their ring membership
  -n             all non ring bonds are marked as "nrings 0"
  -a             only aromatic atoms will match aromatic atoms
  -d             the saturated/unsaturated status of atoms will be preserved
  -s             only allow substitutions at     isotopically labelled atoms
  -w             only allow substitutions at NON isotopically labelled atoms
  -t             not all isotopic atoms need be substituted
  -c             the isotopic number is the number of extra connections at that atom
  -k             use preference values to resolve symmetric atoms
  -u <smarts>    smarts to specify embedding points
  -f ele=smarts  atoms with element type <ele> should match only
                 atoms matching <smarts>
  -h             condense explicit hydrogens to hcount directives on their anchor atoms
  -R <rx>        atoms of type <rx> specify substitution points
                 <rx> is a regular expression, e.g. '^R[0-9]*$', or just 'R'
  -o             remove chirality information from molecules
  -L <smarts>    specify atoms that bind to external group
  -l <nbonds>    include all atoms within <nbonds> of the -L atom(s)
  -I             only include isotopically labelled atoms in the query
  -e             query file to contain just element type and connectivity info
  -V <file>      file containing environment specification
  -X <file>      file containing environment_no_match specification
  -F <fname>     create a file containing the names of all the query files
  -S <fname>     specify output file name
  -x <iso>       atoms with isotope <iso> are translated to match any atom type
  -P <fname>     serialized proto output (use for many queries), 'tsubstructure -q TFPROTO:file ...'
  -p             write individual textproto files
  -b             put all queries in a single file rather than separate file for each
  -D ...         create proto query files with GeometricConstraints
  -Y ...         more obscure options, enter '-Y help' for info
  -i <type>      specify input file type
  -A <qualifier> Aromaticity, enter "-A help" for options
  -g ...         chemical standardisation options, enter '-g help' for info
  -v             verbose operation
```

Fundamental to `mol2qry` is how the information in the molecule gets
converted to substructure query form.

The default behaviour is that the molecule is converted to a form
that can only match itslf. This is likely not useful - but is helpful for
testing. Options
can be specified that relax the constraints put upon the atoms
in the query, allowing the generated query to match more molecules.

Output can be in various forms. The default is as query files, but various
kinds of proto output forms are also supported. The idea is that queries
are generated here, and can then be taken to `tsubstructure` and searches
performed.

## Options

### -m
When a query is constructed from a molecule, each atom has a certain
number of connections and bonds. By default, if an atom starts with
two connections, the query will contain
```
  ncon: 2
```
which means that query atom only matches atoms with exactly two connections.
With the `-m` option, the query will contain
```
  min_ncon: 2
```
which allows atoms with 2 or more connections to be matched.

### -j
Atoms preserve their ring membership. A ring atom that is in 1 ring
will only match an atom that is also in 1 ring. 

### -n
All non ring atoms will have their
```
  nrings: 0
```
atribute set. So, non ring atoms will only match non ring atoms.

### -a
Only aromatic atoms match aromatic atoms. By default, non aromatic
atoms do not have an aromaticity attribute set. But if the `-a` option
is specified, non aromatic atoms only match non aromatic atoms. Without
this, a terminal methyl atom could potentially match an atom in a
benzene ring - assuming that the `ncon` attribute had been relaxed.

### -d
Preserve the saturation status of atoms. This will prevent unsaturated
atoms matching target atoms that are unsaturated.

### -s
One of the most useful options. Only allow substitutions at atoms that
have an isotope. Can be any number. All atoms have their `ncon` attribute
set to what is in the starting molecule. Any atoms with an isotope will
instead have
```
  min_ncon: 
```
specified. That way, substitution can happen only at these atoms.

### -t
By default, there must be substituents at all isotopic atoms. The `-t` option
relaxes that requirement. This is mostly useful when there are multiple
sites in a molecule.

### -c
In this case, the numeric value of the isotope is the number of extra
connections at the atom. So, isotope 1 means that an atom with 2 connections
will only match an atom with 3 connections, but not an atom with 4 connections.

### -I
Only include isotopically labelled atoms in what is converted to query
atoms.

### -k (obscure)
Use preference values to resolve symmetric atoms. If there are symmetric atoms
present, the query generated will not differentiate these atoms. With the `-k`
option, consistency of matches can be gained across multiple matches.

### -u \<smarts\>
Rather than using an isotopic label to specify the atoms where substitutions
are allowed, use a smarts.

### -h
If explicit hydrogens are present, do NOT include those atoms in the
smarts, but set the `hcount` directive for the atoms to which they are
attached. Note that all explicit hydrogens are summed, so `H-C(-H)(-H)` gets
converted to 
```
  min_hcount: 3
```
In retrospect, perhaps it should be the `hcount` rather than `min_hcount`
attribute that gets set. To be investigated.

### -e
Only include element type and connectivity in the query. All bonds become
single bonds. Starting with a benzene will match a cyclohexane.

### -f ele=smarts
Each element of type `ele` in the input molecule can only match `smarts` in
the query. So, if Silicon atoms were being used as an anchor atom, specifying
`-f 'Si=N'` would generate a query where all Silicon atoms only matched a nitrogen.
A more complex use might be `-f 'Si=[ND2x2]'` where Silicon atoms are converted
into query atoms that match a two connected Nitrogen with two ring bonds
attached. Note that only atomic smarts are possible, so changing to matching
a fragment is not possible here. But see ...

### -R \<rx\>
Regular expression match for the element type where attachments are
possible. Often attachment points might be specified by strange elements
like `R1` or `R#`.

### -o
Remove chirality from the input molecules before the query is formed.

### -x \<iso\>
Atoms with isotope \<iso\> are translated to match any atom type.

### -L, -l
These allow creating a subset of the atoms in the incoming molecule
that are converted to query atoms. For example if you were interested
in the different kinds of acids, you might do something like
```
-L '[CD3](-[OH])=O' -l 3
```
which creates a query in which only the atoms that are within 3 bonds
of the first query atom, `[CD3]` are included in the query. Note that
this only identifies the atoms to be included, if you want them to be
able to match more generally than in the query, use the `-m` option
to enable them to match more general connectivity. The first query
atom is matched precisely.

### -V, -X (obscure)
This is used to add environment specifications to the query (-V) or
environment no match specifications (-X). The input
must be a file containing a valid environment specification, in MSI/C2
query form. This is equivalent to forming the query in MSI/C2 format
and then inserting the contents of the file into that query - although
internally the file is interpreted.

At this stage, I do not recall the use case for this. I remember it
involved a complex reaction enumeration where there were important
constraints on what could be used.

## Output
The default output is to create a separate query file, in MSI/C2 format,
for each molecule in the input. By default, those files are named 
`mol2qry0.qry`, `mol2qry1.qry`, ... The file name stem can be changed
via the `-S` option.

### -S \<stem\>
Individual MSI/C2 filec created with names `stem0.qry`, `stem`.qry`,...

### -b
Rather than separate query files, create a single query file containing
all queries.

### -P \<fname\>
Write serialized proto output. This is a binary TFDataRecord file that
can be processed with tsubstructure with `tsubstructure -q TFPROTO:fname ...`

### -p
Write individual textproto query files, subject to the `-S` value.

### -D ...
Create proto query files with geometric constraints. This only works
if the input molecules contain 3D coordinates. The resulting smarts
consists only of the atomic numbers of the atoms, fully disconnected.
So, `CCC` generates a smarts of `C.C.C` which could be very inefficient.
As the need arises, this should be made more flexible.

For each pair of atoms in the molecule, we find the distance between
those atoms, and add a GeometricConstraint specifying that the
distance between the two atoms must be the computed distance, plus or
minus a value specified by `-D tol=...`.

The resulting query should be able to find matches where the conformation
of the matched molecule matches what was in the query. Again, the atom
matching should be made flexible.

## -Y more options.
The `-Y` options are
```
 -Y minextra=n  for a match, target must have at least N extra atoms
 -Y maxextra=n  for a match, target must have at most  N extra atoms
 -Y APPC=<s>    append <s> to the comment field of all queries produced
 -Y exph        add explicit hydrogens, but construct query so anything matched
 -Y ablk        aromatic bonds lose their kekule identity
 -Y minfm=<f>   set the min fraction atoms matched to <f>
 -Y maxfm=<f>   set the max fraction atoms matched to <f>
 -Y A2A=<f>     set aromatic atom translation
 -Y A2A=1       aromatic atoms become 'aromatic'
 -Y A2A=2       aromatic heteroatoms must match aromatic heteroatoms
 -Y A2A=3       aromatic rings must preserve the number of heteroatoms
 -Y rmiso       remove all isotope information from input molecules
 -Y ncon=n      matches must have exactly  <n> connections to unmatched atoms
 -Y min_ncon=n  matches must have at least <n> connections to unmatched atoms
 -Y max_ncon=n  matches must have at most  <n> connections to unmatched atoms
 -Y test        for each query formed, do a match against the starting molecule
```

### minextra=n, maxextra=n
Sometimes it is desirable to avoid having a query match a molecule that is
significantly different in size from the starting molecule. So if you began
with a molecule that contained 6 heavy atoms, and wanted to limit matches to
molecules that contained at most 10 heavy atoms, use `-Y maxextra=4`. The
`minextra` specifies a floor on the number of extra atoms, specifically to
avoid finding what might be uninteresting matches.

### minfm=\<f\>, maxfm=\<f\>
Sets the `min_fraction_atoms_matched` and `max_fraction_atoms_matched`
attributes in the SingleSubstructureQuery proto. This is another way
of controlling how different the molecules matched are from the starting
molecule.

### -Y APPC=\<s\>
Append \<s\> to the comment field of all queries produced.

### -Y exph
Explicit hydrogens are added to the query as query atoms, but they are
instantiated as atoms that can match anything.

### -Y ablk
Aromatic bonds lose their kekule identity. I think this is obsolete now.

### -Y A2A=
Controls how aromatic atoms get translated to query atoms. Possible values
are
1. aromatic atoms become `aromatic: true`
2. aromatic heteroatoms only match aromatic heteroatoms.
3. aromatic rings must preserve the number of heteroatoms
This allows flexible matching of aromatic rings.

### -Y rmiso
Isotopes are removed from the input molecules. Note that this happens
just after reading the molecule, so it will invalidate any attempt at
using isotopes to influence the query formed.

### -Y ncon=n, min_ncon=n, max_ncon=n
Sets the `ncon` attribute in the SingleSubstructureQuery proto. This is
the number of unmatched atoms attached to the matched atoms. This controls
how connected the matched atoms are within the larger molecule

### -Y test
For each query formed, do a match against the starting molecule. This can
be a useful testing tool. Note that it is possible to generate queries that
do not match the starting molecule.
