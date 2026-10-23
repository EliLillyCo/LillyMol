# gfp_group_spread
This tool can be used for selecting groups of molecules.
The anticipated use case is selecting molecules that have been plated, and
selections must be plate-at-a-time.

At each stage of the selection process, the most desirable group is selected.

The desirability of a group can be a function of many different factors.

- The internal diversity of the group.
- The minimum distance between group members and previously selected groups.
- The average distance between group members and previously selected groups.
- The desirability of the molecules in the group.

Every time a group is selected, the distances associated with the remaining
groups are adjusted before the next most desirable group is selected.

For each molecule in the input .gfp file a group membership must be specified. This
can be either via a descriptor file mapping id->group or specifying a column in the
fingerprint name that is the group for that molecule.

## TLDR
```
gfp_group_spread -G file.grp file.gfp
```

where `file.grp` might look like
```
ID Group
ID1 G1
ID2 G1
ID3 G2
ID4 G3
```
and `file.gfp` must contain smiles with names 'ID1', 'ID2', 'ID3' and 'ID4'.

Output might look like
```
G2 0.5 0 ID3
G1 0.3 0 ID1 ID2
G3 0.2 0 ID4
```
where the group selected is in column 1. Then follows a distance to previously selected
measure, an average desirability for that group, and then the group members.

## Options.

### -G group membership

The -G option must specify group memberships. For every fingerprint there must be
a group membership specified. Note that the -G file has a header record, and every
molecule in the fingerprint file must have a group membership specified.

### -D Desirability
Each molecule can have a desirability. Since the average desirability of the molecules in
a group gets combined with things like average distances, it is best for these
measures to also be in the range [0,1], although this is not enforced. The -D file must be a
descriptor file
```
ID Desirability
ID1 0.1
ID2 0.99
ID3 0.80
ID4 0.3
```
Any molecule withouth a desirability value will retain the default value of zero.

### -P Previous Distances
Many times there will be a set of molecules already selected, and so for every molecule
there is a distance to previously selected value available. This can be specified in one
of two ways.

- Supplying a descriptor file of id->distance mappings
- Supplying a .gfp file and distances will be computed.

Beware supplying a large .gfp file, that could result in a lengthy computation.
Note that the output of
```
nplotnn -X table1 file.nn
```
is a suitable input for the -P option, since it is a listing of identifier vs
shortest distance values.

### -w weights
This is the most complex part of this tool.

The desirability of a group of molecules is a function of many factors, as noted above.
Those various factors are set up so that larger values are always more desirable. Each
component can be assigned a relative weight.

If you are going to specify individual weight values, you probably want to start with
```
-w zero
```
which will clear any default weights in effect - a default of 0.5 for average distance
and average desirability.

The following values are recognised.

- avedist=      relative weight for average group distance.
- mindist=      relative weight for shortest distance of any molecule in the group.
- desr=         relative weight for the average desirability of group members.
- intra=        relative weight for the average intra-group distance.

### -T \<dist\>
For standard gfp fingerprints long distances become less meaningful and beyond
a point, two molecules should just be considered "different". The -T option imposes
a cutoff on distance computations. Any computed distance larger than this value is
truncated. Suggest a value like 0.4 or similar. This is not required and the default
is no cutoff.

## Output
The output consists of

- The Group Name
- The average distance to previously selected for that group.
- The average desirability for that group.
- The group members

One of the hardest aspects of this tool is that there is no obviously
correct way of ranking groups at each step. [gfp_spread](gfp_spread.md) uses
a desirability based scaling factor to scale the distance associated with
each candidate fingerprint. This tool attempts to generalise that concept
to the idea of a "candidate" being a set of fingerprints, rather than just one.
