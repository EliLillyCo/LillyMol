# molecule_subset
## Purpose
`molecule_subset` creates a new molecule that consists of just a subset
of the atoms in a starting molecule. The atoms that are copied into the
newly created subset molecule are defined by substructure queries.

This is useful when looking at sets of reagents. Across a set of
reagents, you might ask how many different types of reagent are
there within 3 bonds of the reactive atom?

See also `get_substituents`, a special purpose tool designed to identify
molecule subsets that are defined as all atoms down a bond.


### 3D
`molecule_subset` can operate on 3D molecules, preserving the coordinates.

The `-T` option can include spatially adjacent atoms that are within
a given radius of a matched atom. So, if you wanted atoms that are within
5 Angstroms of an acid, that might be
```
-s [OH]-C=O -T 5.0
```
Note that all matched atoms are considered when determining which
other atoms to include.

The `-C` option can be used to restrict the matched atoms to a region
of space - if you have aligned molecules. Specify one or more spatial
restrictions with the -C option, and then only matched atoms that
satisfy these constraints will be used for subsequent calculations.
The previous example of the acid match, could be restricted to only
those acid groups in a region of space

```
-s [OH]-C=O -T 5.0 -C 5.0,6.5,-3.2,8
```
which means that all matched atoms must be within 8.0 angstroms of
the point (5.0, 6.5, -3.2). Note that it is to be expected that
this may result in all matched atoms no longer matching.

## HOWTO
The usage message is
```
Creates subsets of molecules defined by query matches
  -q <query>     specify substructure query
  -s <smarts>    specify smarts (ss)
  -n             write the molecule subset which does not match the query
  -b             create the subsets by bond rather than by atom
  -m <number>    only process the first <number> hits for any query
  -m do=<nn>     process embedding number <nn> (starts with 0)
  -m each        process each embedding separately
  -z i           ignore molecules not matching any query
  -z w           write molecules not matching any query
  -z b           blank line to stdout for non-matching molecules
  -z nmapp=xxx   append 'xxx' to the names of non-matching molecules
  -f             also write small fragments with selected atoms
  -T <dist>      with 3D input, include atoms within <dist> of any matched atom
  -C x,y,z,rad   with 3D input, only search atoms within <rad> of <x,y,z>
  -i <type>      specify input file type
  -o <type>      specify output file type(s)
  -S <string>    create output files with name stem <string>
  -M ...         various miscellaneous options, enter '-M help' for info
  -E ...         element specifications, enter '-E help' for details
  -A ...         standard aromaticity options, enter '-A help' for details
  -v             verbose output
```
Most of the options are unsurprising.

Using composite queries, `C&&N` is tricky, because the matched atoms
reported will be those atoms that were matched by the last of the
individual queries to match. So `C&&N` will report the Nitrogen atom,
whereas `C||N` will report the Carbon (TODO:ianwatson seems to not
be the case, investigate).

The `-b` option is interesting. By default, all matched atoms are included
in the subset. If the `-b` option is used, the selection is based on
bonds rather than atoms. Within the matched atoms, all atoms that are
bonded will be included in the subset. So, if there is a disconnected
atom in the query, that will not be included. To find acids that
have a fluorine atom < 10 bonds away
```
-s '[OD1]-C(=O)...{<10}F' -b
```
should work. The subsets generated are just the acid atoms.
