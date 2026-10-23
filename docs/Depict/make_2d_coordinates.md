# make_2d_coordinates

Generate 2D coordinates for depiction.

LillyMol reads and writes 3D coordinates but has never generated 2D layouts of
its own, so molecules coming from smiles had nothing sensible to draw. This tool
computes a layout — the flat, even-bond-length drawing a chemist expects — and
writes it as a molfile.

The layout engine is Schrodinger's
[coordgen](https://github.com/schrodinger/coordgenlibs), the same library RDKit
uses for its high quality depictions. The library API is
`depict::Generate2DCoordinates` in `src/Depict/coords_2d.h`.

## Building

coordgen is an optional dependency, so this tool is **not built by default**. Set
`BUILD_SCHRODINGER_2D` before building and `src/build_linux.sh` will fetch and
build coordgen into `third_party/`, then build and test the `Depict` package
against it:

```
export BUILD_SCHRODINGER_2D=1
make
```

Everything else in LillyMol builds and runs exactly as before without it.

## Usage

```
make_2d_coordinates -S out file.smi
```

Only `.sdf` can hold coordinates, so that is the default output type — there is
no need to give `-o sdf`.

```
 -b <length>    bond length in the generated layout, default 1.5, the molfile convention
 -p quick       lay out quickly, accepting clashes that would otherwise be resolved
 -p standard    the default effort
 -p best        work hardest to resolve clashes, slowest
 -e             spread substituents evenly around an atom
 -k             skip the force field refinement that follows construction
 -n             do not centre the layout on the origin
 -t             write LillyMol's terse molfile counts line instead of full V2000
 -w             do not assign wedge bonds
 -h             remove explicit hydrogens first
 -c             discard molecules for which no layout could be generated,
                rather than treating that as a fatal error
 -l             reduce to the largest fragment
 -g ...         chemical standardisation options
 -S <fname>     output file name stem
 -i ...         input file specification
 -o ...         output file specification. Default .sdf
 -v             verbose output
```

## Full V2000 records

Molfiles are written with **full V2000 counts lines** rather than LillyMol's
default abbreviated form, which some toolkits refuse to read. Coordinates only
LillyMol can read would defeat the purpose, so this is the default here even
though it differs from the rest of LillyMol. Use `-t` for the terse form.

## Explicit hydrogens

Explicit hydrogens are laid out as atoms in their own right, which is usually
not what a depiction wants — the drawing ends up dominated by hydrogens. Use
`-h` to remove them first. This happens after standardisation, so `-g` options
that add or remove hydrogens are respected.

```
make_2d_coordinates -h -g all -S out file.smi
```

## Choosing an effort level

`-p` controls how hard coordgen works to separate atoms that would otherwise
overlap. For ordinary drug-like molecules all three settings usually give the
same answer, so the default is a reasonable choice; `-p best` is worth trying on
crowded fused or bridged ring systems, and `-p quick` on very large inputs where
throughput matters more than the last few overlaps.

As a rough guide to cost, a 39,823 molecule file took about 4 minutes at the
default precision on a single core.

## Failures

A molecule that cannot be laid out is a **fatal error** by default, since it is
rare enough that it usually means something unexpected. `-c` skips such
molecules instead and reports a count under `-v`.

Failures are of three kinds:

- an empty molecule;
- a single connected fragment containing more than 40 rings, a coordgen limit;
- coordgen produced NaN coordinates, or collapsed part of the structure onto a
  single point. Both are only seen on heavily bridged polycyclics.

The last case matters because coordgen reports these coordinates as successfully
set. Writing them out would produce a molfile no toolkit can read, or one with
atoms drawn on top of each other, so they are detected and reported as failures
instead — a molecule that is written is always drawable. Over a 39,823 molecule
structure-mutation set 47 molecules failed this way; over 1,893 ChEMBL
molecules, none did.

## Layout quality

A successful layout is guaranteed to be *drawable*, not necessarily *pretty*.
Bond lengths are uniform for ordinary molecules but cannot be for bridged ring
systems, where the bridge bonds have to be drawn short — that is the standard
2D convention, not a defect. Over 1,893 ChEMBL molecules, 89.9% of bonds came
out within 5% of the requested length, and 99.8% of molecules kept every
non-bonded pair of atoms at least half a bond length apart.

With `-v` the tool reports how many layouts have non-bonded atoms closer than
half a bond length. These are still written; they will simply look crowded.

## Stereochemistry

Chiral centres are given **wedge and hash bonds**, so chirality survives being
handed to other software. This matters because LillyMol also writes chirality in
the MDL atom parity field, and while it reads that back itself, other toolkits
ignore parity on a 2D record: without wedges a depiction is stereochemically
silent no matter how correct the molecule is inside LillyMol. RDKit reads back
every stereocentre in a test set covering acyclic, ring, quaternary, bridged and
sugar centres.

Wedges are assigned after the layout, since the direction one points depends on
the coordinates, and any wedge the input already carried is discarded — it was
drawn against different coordinates and is no longer a true statement. `-w`
turns the whole thing off, leaving parity as the only record of chirality.

Under `-v` the chiral centres are reported in three groups:

```
25 chiral centres, 20 given a wedge bond
5 marked chiral but not stereogenic, nothing to draw
```

The middle group is the normal case of an atom marked chiral in a smiles that is
not in fact stereogenic — `C[C@H]1CC[C@@H](C)CC1`, whose two ring branches from
each centre are identical. There is no stereochemistry to draw, so none is
drawn. A third line, *"stereocentres could not be wedged"*, is the case that
actually loses information; it does not appear for ordinary structures.

Cis/trans (E/Z) is **not** depicted yet. It is preserved — the single bonds that
carry it are never taken for a wedge — but nothing is written that conveys it to
a reader working from the drawing.

The library API is `depict::AssignWedgeBonds` in `src/Depict/wedge_2d.h`. It
needs only coordinates that are already present, so it can be used on molecules
laid out by something else, and it does not pull in coordgen.

## See also

- `src/Depict/coords_2d.h` — the library interface, for calling the layout
  directly rather than through this tool.
- `src/Depict/wedge_2d.h` — wedge bond assignment on its own.
- `src/build_linux.sh` — the `BUILD_SCHRODINGER_2D` block that builds coordgen.
