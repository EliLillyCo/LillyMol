# Grid Proximity

This tool takes a grid produced by Insight, and one or more
ligands. For each ligand atom, it determines which grid points
are within a fixed radius of that ligand atom. Once the grid
points are identified, the result is the sum of the data values
present at each grid point.

## Options

The usage message is
```
 -G <fname>     file containing InsightII grid
 -t <dist>      count density within <dist> of matched atoms in the ligands
 -q <qry>       query specifying ligand atoms to score
 -s <smarts>    smarts specifying ligand atoms to score
 -z i           ignore molecules not matching the query
 -M <fname>     write the grid as a molecule
 -m <float>     multiple grid raw values by <float> and place isotopes
 -i ...         input specification(s)
 -A ...         aromaticity options
 -E ...         element related options
```

If only a subset of the atoms in the molecule are to be processed,
that can be specified via the `-q` or `-s` options. If such queries
are specified, then the `-z i` option tells it to ignore any molecules
that do not match those queries, by default it will die.

The `-M` option is a useful diagnostic. It will write the ligand, and
all the grid points except those that are close to the ligand. When
viewed you should see the ligand, surrounded by an empty shell, and
then grid atoms.

The properties on the grid are floating point numbers. It can be useful
to visualise the grid properties as isotopic labels. But if those grid
properties are small floating point numbers, it may be necessary to
multiply them by a factor in order for them to show up as useful isotopic
numbers. Use the `-m` option to do this.

Again, the -M and -m options are mostly useful for diagnostic work.

Today the tool works as a fixed distance. It would not be hard to make
it use the vdw radius - or any radius that was wanted.

And if the grid were weighted, with the most important parts of space
given higher weights, then a meaningful score could be derived.

Execution times were trivial on the cases I tried.

## Extensions
The tool as it exists just compares molecules individually to a fixed
grid. It does not compare ligands themselves.
 
It would be easy to alter the tool to generate a fingerprint
where a bit would be set if the ligand is close to that grid point. In
fact I just did that, but I do not have an Insight grid file or ligand
to test.

I have set it up to do both binary fingerprints and sparse fingerprints.
Sparse fingerprints could also encode the grid property (weight), but
that is not in place now - but trivial to do. But I suspect that
using sparse fingerprints in calculations might be slow, although I
don't really know...

But what this should do is generate a .gfp file entry for each ligand
where a bit is set if it is close to that ligand. Then all the usual
similarity tools become available. Depending on the resolution of the
grid, this might be quite large, so hard to know what to expect.
