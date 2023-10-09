# rgroup

`rgroup` is one of a series that identifies substituents at specific points
in a molecle, or in a series of molecules. It can be the precursor for an R group
analysis. See also `remove_and_label`, `get_substituents` and others.

The tool is complex since it evolved over time to satisfy several different needs.

## Usage
The following options are recognised
```
Must specify one or more queries via -q or -s options
Molecule_Tools/rgroup.cc compiled 2023-Jun-20 git hash 944f447
Usage: ./bazel-bin/Molecule_Tools/rgroup -q/-s query <options> <input file>
Identifies substituents based on what is attached to query atoms
  -q <file>      specify query(s)
  -s <smarts>    queries as smarts
  -r q.n         find R group from atom N in query Q
  -z i           ignore molecules which don't match a query
  -u             unique embeddings only
  -d <ele>       add an embedding dummy atom, type <ele>, to all groups
  -a             process all embeddings
  -b             process symmetry equivalent embeddings
  -e             do multiple substituents on same atom separately
  -k             include attachment atom in R groups
  -I <iso>       label join points with isotope <iso>
  -X <symbol>    after building, remove all elements of type <symbol>
  -c <natoms>    min atoms in an R group
  -C <natoms>    max atoms in an R group
  -H             make implicit hydrogens explicit
  -h             after processing, remove explicit Hydrogens not attached to isotope (-I)
  -y             process substituents involved in ring bonds
  -i <type>      input type
  -A <qualifier> Aromaticity, enter "-A help" for options
  -K <...>       Enter "-K help" for SMILES options
  -E <symbol>    create element with symbol <symbol>
  -v             verbose output
```

Fundamental to the usage is specification of the query, or queries, that define how the
substituents are attached to the scaffold. Queries are defined via the `-s` option (smarts)
or `-q` options, any of the usual `-q` option qualifiers.

In addition, the actual substitution points must be specified via the -r option. Substitution
points are defined by two numbers, the query number (q) and the matched atom number (n). Both
of these start at zero. So `-r 0.2` would look for substituent groups that branch off from
matched atom 2 in query 0. Any number of -r options can be specified.

Note that if no `-r` option is specified, all matched atoms are used as starting
points for a substituent.

If multiple queries are used, they must all be consistent - they must have the
same number of query atoms. It is assumed that what is found as a substituent
at each matched atom is the same across queries, although that is not enforced.
You may need to use a query environment, or recursive smarts, in order to get
each query to have the same number of atoms, while aligned.

### Query Atom Alignment.
For example we are looking for substituents attached to either an amide, or a sulfonamide.
That would be done with
```
-s 'O=[$([CD3]),$([SD4]=O)]-[ND>1]' -a
```
or
```
-s 'O=[CD3]-[ND>1]' -s 'O=[$([SD4]=O)]-[ND>1]' -a
```
yield the same structures - although different naming schemes.

A query file for the sulfonamide can use the query environment to specify atoms
that must be matched, but which will not be part of the embedding.
```
query {
  smarts: "O=S-[ND>1]"
  environment {
    attachment {
      attachment_point: 1
      btype: SS_DOUBLE_BOND
    }
    smarts: "[OD1]"
  }
}
```
Again, the alignment with the regular Amide query is preserved. This time the `O=`
group will be randomly assigned.

If one of the queries does not match, the tool fails, unless `-z i` is specified,
in which case it ignores molecules not matching.

## -a process all embeddings
By default, only the first embedding of the substructure in the molecule is processed.
With this option, all matches are processed. The resulting substituent names will reflect
differing embedding numbers.

## Matching Control.
The `-u`, unique embeddings only, and -b (process symetry related matches) do what
theyusually do. See (tsubstructure)[tsubstructure.md].

## Labelling
You almost certainly want to label the attachment points with the `-I` option.

## Size
the `-c` and `-C` options control the size of the substituents that are reported.

## Hydrogens
If the `-H` option is specified, all implicit Hydrogens are converted to explicit. After
processing, these can be removed with the `-h` option, except that any Hydrogen that is
attached to an isotopically labeled atom (an attachment point) will not be removed.


## Multiple Substituents
There is ambiguity as to what should happen if there are multiple substituents
at an anchor atom. For example in the (very contrived) example of `C1CCC1(C)O`
with the query
```
-s '[CD4]1CCC1'
```
generates output
```
C1CCC1(C)O p2
O[1CH2]C RG 0.0.0.0
```
where the anchor atom is included with the substituent. With the `-e` option,
which processes each substituent separately, this becomes
```
C1CCC1(C)O p2
[1CH4] RG 0.0.0.2
[1OH2] RG 0.0.0.3
```
where now there are two separate substituents.

## Output Format
The output is a smiles file that might look like (the query was `c1ccccc1`)
```
C1(=C(C=CC(=C1O)O)[C@H](O)CN)F CHEMBL40796
[1FH] RG 0.0.0.0
O[1CH2]CN RG 0.0.1.0
[1OH2] RG 0.0.4.0
[1OH2] RG 0.0.5.0
```
All substituents have an 'RG' as the first token of the name.

The sequence of digits is derived from
```
query_number.embedding_number.query_atom.index
```
In the example above there was only one query, so all the first digits are zero.

Again, since the `-e` option was *not* used, there is only one embedding, so 
the second digit is also all 0.

Substituents were found at three positions on the ring, that their indices show
up as the third digit.

And since there were no multiple matches at those points, the last digit is
also zero.

## Other options.
### -d \<element\>
Rather than placing an isotope at the attachment points, it make make sense to
add a new atom instead.

### -y
Generally the tool avoids substituents that involve ring bonds at the
attachment point.

For example, by default the molecule `CCC1CC1` when processed with
```
rgroup.sh -s '[CD1]CC' -r 0.2 ...
```
generates no output. It detects that there are ring bonds on the attachment
point. However if the `-y` option is used, the following output is generated
```
CCC1CC1 t1
C1CC1 RG 0.0.0.0
```
Care is needed with this option. With the `-y` option in effect. If the matched
atom is in a ring, those paths will not be followed. Perhaps this should
be the default.

