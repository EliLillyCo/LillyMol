# Molecular Grid

Molecular grid creates a grid of points that encompas a set
of ligands. The usage message is
```
  -d <dist>     spacing between grid atoms
  -x <dist>     ensure grid expands <dist> from any atom (default 2A)
  -t            trim the grid to only those atoms with -x of a ligand atom
  -j            for each ligand, compute total and min grid point occupancy
  -r <number>   during the -t option, report progress every <number> ligands examined
  -S <fname>    write grid to <fname>
  -o <type>     output type (default sdf)
  -e <ele>      element to use for the grid
  -l            reduce to largest fragment
  -i <type>     input specification
  -g ...        chemical standardisation options
  -E ...        standard element specifications
  -A ...        standard aromaticity specifications
  -v            verbose output
```
The primary parameters governing the operation are the `-d` and `-x`
options.
