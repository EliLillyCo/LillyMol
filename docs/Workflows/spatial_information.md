# Spatial Information.

Most LillyMol tools are 3d capable. There are some specific tasks for which
specific 3D tools have been built. As is often the case with tools developed
over time, there can be overlapping functionality.

## get_coordinates.
Given a query specification, smarts or query file, this tool extrracts the coordinates of matched atoms.
A variety of output forms are available. For example to extract the coordinates of a phenolic
oxygen atom one might do something like
```
get_coordinates -s '[OH]c1ccccc1' file.sdf
```
Might produce output such as
```
C1(=C([H])C(=C([H])C(=C1[H])[H])[H])O[H] CHEMBL14060 0 0
0 1 2.2673 0.0194 0.0148 O [OD2H1v2]
1 0 0.904 -0.0054 0.0111 C [cD3H0v4;r6]
2 2 0.2004 -1.2066 0.0063 C [cD3H1v4;r6]
3 3 -1.1955 -1.1829 -0.0026 C [cD3H1v4;r6]
4 4 -1.8752 0.0364 -0.0068 C [cD3H1v4;r6]
5 5 -1.1601 1.2349 0.0008 C [cD3H1v4;r6]
6 6 0.2344 1.2132 0.0103 C [cD3H1v4;r6]
-----------------
```
Note that the coordinates of all matched atoms are included in the output.

The first column is a sequential index corrdsponding to the matched
atom. This enables comparability across molecules. 

The second column is the atom number within the original molecule, which may
or may not be of interest.

Then follows the coordinates of those atoms, the atomic symbol and a
possible smarts for that atom.

Note that in this case the query is symmetric, so in addition to the output
above we also see
```
C1(=C([H])C(=C([H])C(=C1[H])[H])[H])O[H] CHEMBL14060 0 1
0 1 2.2673 0.0194 0.0148 O [OD2H1v2]
1 0 0.904 -0.0054 0.0111 C [cD3H0v4;r6]
2 6 0.2344 1.2132 0.0103 C [cD3H1v4;r6]
3 5 -1.1601 1.2349 0.0008 C [cD3H1v4;r6]
4 4 -1.8752 0.0364 -0.0068 C [cD3H1v4;r6]
5 3 -1.1955 -1.1829 -0.0026 C [cD3H1v4;r6]
6 2 0.2004 -1.2066 0.0063 C [cD3H1v4;r6]
-----------------
```
which shows the matched atoms in a different order. If the only atom of interest
was the Oxygen, these other matches would not be of interest, and the `-u` option
could be used to suppress matches which are not unique in 2D. But note that matches
that are symmetric in 2D will be distinct in 3D.

It may also be of interest to add a centroid atom defined by all matched atoms, `-C Ce`
```
C1(=C([H])C(=C([H])C(=C1[H])[H])[H])O[H] CHEMBL14060 0 0
0 1 2.2673 0.0194 0.0148 O [OD2H1v2]
1 0 0.904 -0.0054 0.0111 C [cD3H0v4;r6]
2 2 0.2004 -1.2066 0.0063 C [cD3H1v4;r6]
3 3 -1.1955 -1.1829 -0.0026 C [cD3H1v4;r6]
4 4 -1.8752 0.0364 -0.0068 C [cD3H1v4;r6]
5 5 -1.1601 1.2349 0.0008 C [cD3H1v4;r6]
6 6 0.2344 1.2132 0.0103 C [cD3H1v4;r6]
. . -0.08924288 0.01557144 0.004842857 Ce .
-----------------
```
The centroid will be the same for all matches.

Using `-o csv` will result in csv output which might be easier process. It is also
possible to obtain textproto output.

## molecule_subset
Rather than writing coordinates, this tool writes a complete molecule containing
just the atoms specified by the queries. For example using the phenol example from
above
```
molecule_subset -m each -M isomatch -z i -s '[OH]c1ccccc1' -S /tmp/h -o sdf file.sdf
```
will generate outputs that might look like
```
CHEMBL14060
Blank
Blank
  7  7
    0.9040   -0.0054    0.0111 C   0  0
    2.2673    0.0194    0.0148 O   0  0
    0.2004   -1.2066    0.0063 C  -9  0
   -1.1955   -1.1829   -0.0026 C  -8  0
   -1.8752    0.0364   -0.0068 C  -7  0
   -1.1601    1.2349    0.0008 C  -6  0
    0.2344    1.2132    0.0103 C  -5  0
  1  2  1
  1  3  2
  1  7  1
  3  4  1
  4  5  2
  5  6  1
  6  7  2
M  ISO  7   1   2   2   1   3   3   4   4   5   5   6   6   7   7
$$$$
```
Here the atoms are not written in the order of the query, but in the order of the
starting molecule. But because we added `-M isomatch` the atoms are labelled
by the query atom number (+1) of the query atom that matched. Therefore atom 2,
the Oxugen atom, has isotope 1.

molecule_subset is a complex utility, and is more fully described at
[molecule_subset](/docs/Molecule_Tools/molecule_subset.md).

## superimpose_by_matched_atoms
This tool performs alignment to an exising 3D template molecule. Usually this is
used when a known pose is to be applied to a series of unaligned molecules. Note
that the superimposition is rigid. There are other tools that can do flexible
alignments.

The aligned template molecule is specified via the -T option, and the matched
atoms to be considered must be specified as either smarts or query file. The
tool can compute the RMS before superimposition takes place as well as after.

The following invocation demonstrates some of the functionality
```
superimpose_by_matched_atoms -r -R RMS.txt -Z all -s '[OH]c1ccccc1' -Z i -Z ignorezero -S superimposed \
        -o ISIS -z i -z all -v -T template.sdf file.sdf
```
