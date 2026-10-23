# tsubstructure

`tsubstructure` is a command line substructure searching tool. It reads
queries and one or more files of molecules, and reports which molecules match. It can also
write the matching and non matching molecules, label the atoms that matched,
generate fingerprints or descriptors, and
produce per query counts. Most substructure searching in LillyMol starts here.

This page describes the options that are used day to day, with examples that
were run against the current code. The tool has many more options, and the
built in help is complete: run `tsubstructure` with no arguments, or
`tsubstructure -m help`, `-M help`, `-j help`, `-i help`, `-g help`, `-A help`,
`-P help`. How to write queries is described in
[Substructure Searching](../Molecule_Lib/substructure.md), and the equivalent
Python interface in [tsubstructure (python)](../python/tsubstructure.md).

## Contents

* [Quick start](#quick-start)
* [Specifying queries](#specifying-queries)
* [Reading molecules](#reading-molecules)
* [Writing results](#writing-results)
* [What counts as a match](#what-counts-as-a-match)
* [Several queries](#several-queries)
* [Labelling matched atoms](#labelling-matched-atoms)
* [Aromatic bonds and Kekule forms](#aromatic-bonds-and-kekule-forms)
* [Per query results](#per-query-results)
* [Debugging a query](#debugging-a-query)
* [Preparing molecules before searching](#preparing-molecules-before-searching)
* [Behaviour worth knowing](#behaviour-worth-knowing)
* [Option summary](#option-summary)

## Quick start

Count the molecules that contain a nitrile
```
tsubstructure -s 'C#N' file.smi
```
Nothing is written to stdout. The result is reported on stderr
```
9 molecules read, 2 molecules match fraction 0.222222
```
To keep the molecules that matched, use `-m`
```
tsubstructure -s 'C#N' -m nitriles file.smi
```
which writes `nitriles.smi`. To keep the ones that did not match as well, add `-n`
```
tsubstructure -s 'C#N' -m nitriles -n not_nitriles file.smi
```
Use `-m -` to write the matches to stdout, which is convenient in a pipeline
```
tsubstructure -s 'C#N' -m - file.smi | head
```
The input file can be `-`, meaning stdin.
```
cat file.smi | tsubstructure -s 'C#N' -
```

## Specifying queries

A SMARTS is given with `-s`. It can be followed by a name, which is then used in
the per query output described below
```
tsubstructure -s 'C#N nitrile' file.smi
```
The `-s` option can be repeated, and a molecule matches if it matches any of the
queries (see [Several queries](#several-queries)). For more than a few queries, use a file.

| Option | Meaning |
| --- | --- |
| `-s <smarts>` | a SMARTS, optionally followed by a name |
| `-q <file>` | a query file. A historical query file, or a textproto query |
| `-q S:<file>` | a file of SMARTS, one per line, each optionally followed by a name |
| `-q F:<file>` | a file that contains the names of query files |
| `-q M:<file>` | a file of molecules. Each molecule becomes a query |
| `-q PROTO:<file>` | a file containing a textproto query |
| `-q PROTOFILE:<file>` | a file containing the names of textproto query files |
| `-q 'proto:query { smarts: "n" }'` | a textproto query written on the command line |

For example, with `queries.smt` containing
```
C#N nitrile
[SX4](=O)=O sulfonyl
```
the command
```
tsubstructure -q S:queries.smt file.smi
```
matches molecules containing either. Whether a `-q` file is a historical query or a
textproto is usually detected from the file contents, and `PROTO:` removes any doubt.
The textproto form is described in [substructure_proto.md](../Molecule_Lib/substructure_proto.md),
and the query files supplied with LillyMol are in `data/queries`.

## Reading molecules

The input is a smiles file by default, in the form `smiles name`. Other formats
are chosen with `-i`, for example `-i sdf`; `-i help` lists the many input
qualifiers. Two are useful when testing a query on a large file

| Option | Meaning |
| --- | --- |
| `-i do=N` | only process the first N molecules |
| `-i skip=N` | skip the first N molecules |

`-M minat=N` and `-M maxat=N` restrict the search to molecules with at least or
at most N atoms. When both are given, the minimum must be smaller than the maximum.

Molecules are searched as read. See
[Preparing molecules before searching](#preparing-molecules-before-searching) for options that
change them first.

## Writing results

### Which molecules, and in which format

`-m <stem>` writes the molecules that matched and `-n <stem>` those that did not.
The file name is the stem plus a suffix for the output type, which is `.smi` unless
`-o` says otherwise
```
tsubstructure -s 'C#N' -m hits file.smi               # hits.smi
tsubstructure -s 'C#N' -m hits -o sdf file.smi        # hits.sdf
tsubstructure -s 'C#N' -m hits -o smi -o sdf file.smi # hits.smi and hits.sdf
```
If the stem already ends in the suffix for the output type, it is not added a second
time, so `-m hits.smi` writes `hits.smi`. A stem with some other extension has the suffix
added, so `-m hits.txt` writes `hits.txt.smi`, and `-m hits.smi -o sdf` writes
`hits.smi.sdf`. Only a suffix that follows a period is recognised, so `-m out_smi` writes
`out_smi.smi`.

### Extra information on each molecule

Extra `-m` qualifiers add details of the queries that matched. With
`tsubstructure -s 'c q1' -s 'C#N q3' -m - -m QDTVB`
```
C1=CC=CC=C1C1=CC=CC=C1 biphenyl |12 matches q1|
C1=CC=CC=C1 benzene |6 matches q1|
CC#N acetonitrile |1 matches q3|
```
| Qualifier | Adds |
| --- | --- |
| `-m QDT` | `(4 matches to 'c')`, readable |
| `-m QDTVB` | `\|4 matches q1\|`, easier for other programs to parse |
| `-m CSR` | `q1:4 q2:1`, compressed. `-m CSR:,` uses `,` as the separator, and `-o csv` writes a true csv file |
| `-m SEPA` | writes the matches to each query to a separate file, `FOO0.smi`, `FOO1.smi`, ... for `-m FOO` |
| `-m NONM`, `-m NONMX` | also write the non matches to the `-m` file. `NONMX` does not append the query details to them |

The counts are numbers of embeddings, see [What counts as a match](#what-counts-as-a-match).
`-M app=<text>` appends `<text>` to the name of every molecule that matches.

### Other output

| Option | Result |
| --- | --- |
| `-R <file>` | a report of how many molecules each query matched. `-R idfirst -R <file>` puts the query name first |
| `-G <file>` | the atom numbers that matched in each molecule, numbered from zero, `acetonitrile (1 2)` |
| `-J <tag>` | Daylight style fingerprints in TDT form, with `-y <nbits>` setting the size |
| `-a` | a table of per query match counts, see [Per query results](#per-query-results) |

Note that with fingerprint output, if the tag starts with 'FP', the output will be a fixed
width, binary fingerprint. If the tag starts with 'NC' it will be a sparse, counted fingerprint.

## What counts as a match

A molecule matches a query if the query can be placed on it at least once. The numbers of
matches that are reported are numbers of *embeddings*, the distinct ways of placing
the query, and symmetry means there can be more than you expect. With the query `c1ccccc1`
on benzene (`-m - -m QDT`)

| Options | Reported |
| --- | --- |
| none | `(12 matches to 'ring')` - six starting atoms, two directions |
| `-u` | `(1 matches ...)` - unique matches only, ignoring those that cover the same atoms |
| `-f` | `(1 matches ...)` - stop after the first embedding |
| `-k` | `(1 matches ...)` - do not perceive symmetrically equivalent matches |
| `-r` | `(6 matches ...)` - one embedding for each atom that the root atom of the query matches |
| `-M maxe=N` | stop after N embeddings for each query |

This is why a symmetric group such as `[SX4](=O)=O` reports `2 matches` on a sulfone
without `-u`. It does not change which molecules match, only the counts. Numeric
qualifiers in the SMARTS, such as `1[OH]-C=O`, count embeddings in the same way,
see [Substructure Searching](../Molecule_Lib/substructure.md#numeric-qualifiers).

`-M print` prints the atoms of every embedding found (numbered from zero),
which is useful when working out what a query is doing.

## Several queries

By default, a molecule matches if it matches *any* of the queries.
With `-M mmaq`, it must match *all* of them
```
tsubstructure -s 'C#N' -s 'c' file.smi             # a nitrile, or an aromatic atom
tsubstructure -M mmaq -s 'C#N' -s 'c' file.smi     # a nitrile, and an aromatic atom
```
Two options stop work early, and they behave differently

* `-b` stops testing a molecule against the remaining queries after the first one that
  matches. The set of molecules that match is the same as the default, but the details
  reported (`-m QDT` and so on) list only the first query that matched.
* `-B` stops after the first query that does *not* match. A molecule that fails a
  query is not tried against the ones after it. The result therefore depends on the
  order of the queries, and it is not a way of requiring all queries to match. For
  that, use `-M mmaq`. Neither option can be combined with `-M mmaq`.

For instance with `-s 'C#N'` and `-s 'c'` on a biphenyl, an acetonitrile and a phenylacetonitrile,
`-B` matches the two nitriles, and with the queries the other way round it matches the
biphenyl and the phenylacetonitrile.

## Labelling matched atoms

`-j` places isotopic labels on the atoms that matched, which is how
`tsubstructure` is used to mark positions for other tools. The labelled molecules are
written by `-m`. With the queries `-s 'S(=O)=O sulf' -s 'C#N nit'`, shown on
dimethylsulfone and acetonitrile

| Options | dimethylsulfone | acetonitrile |
| --- | --- | --- |
| `-j 1` | `C[1S](=[2O])(=[3O])C` | `C[1C]#[2N]` |
| `-j 1 -j same` | `C[1S](=[1O])(=[1O])C` | `C[1C]#[1N]` |
| `-j qnum` | `C[1S](=[1O])(=[1O])C` | `C[2C]#[2N]` |
| `-j iquery` | `CS(=[1O])(=[2O])C` | `CC#[1N]` |

The number given to `-j` is the first label, and by default the label increases
with each matched atom. The combination `-j 1 -j same`, giving every matched atom the
same label, is the one that is most often wanted. `-j qnum` labels by the number
of the query, counting from one. `-j iquery` labels by the position of the atom
in the query, counted from zero, and an atom with label zero carries no isotope, so
the first query atom is unlabelled.

Molecules that match nothing are not written unless `-m NONMX` is added, which
is usually what is wanted with `-j`.

More forms exist, for example `-j writeach=<file>` writes every labelled embedding
as a separate molecule, `-j amap` uses atom map numbers instead of isotopes, and
`-j noclear` keeps existing isotopes on atoms that did not match. See `-j help`.

## Aromatic bonds and Kekule forms

A single bond in a query matches the single bonds inside aromatic rings, which is
different from the behaviour of most other toolkits. LillyMol keeps the Kekule bond
orders of an aromatic ring alongside its aromaticity, and by default a query bond can
match either. As a result `a-a`, which most people take to mean two aromatic atoms
joined by a single bond outside a ring, matches every molecule that has an aromatic ring.

Run on benzene, naphthalene, pyrrole, furan, 2-pyridone, biphenyl and
2-phenylpyrrole, the query `a-a` gives

| Option | Molecules that match `a-a` |
| --- | --- |
| none, or `-M kekule` | all seven |
| `-M nokekule` | biphenyl, phenylpyrrole |
| `-M fkekule` | all but benzene |

* `-M kekule` is the default. An aromatic bond matches aromatic query bonds, and also
  single and double query bonds through its Kekule form.
* `-M nokekule` makes an aromatic bond match only aromatic query bonds. Then `a-a` is
  two aromatic atoms joined by a bond that is not aromatic, as in RDKit.
* `-M fkekule` is in between: aromatic rings that have alternating Kekule forms, such
  as benzene, lose them, and others keep them. Which rings are treated this way is
  decided by LillyMol and not by the query. For the example above only benzene
  changed.

The behaviour is the same for `c-c`, `a=a` and any other query that puts a single or double
bond between aromatic atoms. The aromatic bond `:` is unaffected.

The simplest way to avoid the difference is to write the query so that it does not depend
on the mode. A bond between two aromatic systems is `a-!@a`, a single bond that is not in
a ring, and it gives the same answer with all three options.

The setting is global, so it applies to all the queries in a run. Queries that you did not write
may depend on the default. On a sample of 200,000 molecules, about thirty of the query files in
`data/queries` give different answers with `-M nokekule`. Some depend on the default deliberately:
`hbonds/imine.qry`, one of the H-bond acceptor queries, finds aromatic ring nitrogens through their
Kekule double bonds, and several of the charge assigner and medchem rule queries are written
in the Kekule form. Others, such as the boronic acid, coumarin, `hydrazide_cyclic` and `imidazole_basic`
queries, have been written so that they give the same answer in every mode. Check the queries that you
rely on before changing the setting. The script `contrib/bin/kekule_query_audit.sh` runs a set of query files
over a set of molecules in both modes and lists the queries that differ, see the
[contrib/bin](/contrib/bin/AAREADME.md#kekule_query_auditsh) notes. In Python this is
`set_aromatic_bonds_lose_kekule_identity()`, see
[tsubstructure (python)](../python/tsubstructure.md#aromatic-bonds-and-kekule-forms).

## Per query results

When there is more than one query, `-a` writes a table with a row for each molecule and a
column for each query, containing the number of matches
```
tsubstructure -s 'c q1' -s 'C#N q2' -a file.smi
```
```
Name q1 q2
biphenyl 12 0
benzene 6 0
sulfonamide 0 0
```
The rows are named by the molecule names in the input file, so without names the first column is
empty. `-Y <stem>` names the columns with the stem and a number, `XX0 XX1`, in place of the
query names. `-M owdmm` writes rows only for molecules that match, and `-M anmatch` writes
a single column, `Matches`, with the total number of matches.

`-v` reports how many molecules matched each number of queries, and the hits for each query
```
3 molecules matched 0 queries
5 molecules matched 1 queries
1 molecules matched 2 queries
0 Details on hits for query 'q1'
Tested 9 molecules, 5 matched, 4 did not. Fraction 0.555556
 4 molecules had 6 hits
 1 molecules had 12 hits
```
The `-R` option writes the numbers of molecules that matched each query to a file.

## Debugging a query

When a query does not match what you expect, first reduce it to something that does, and build
up from there, checking each step on a few molecules chosen so that you know the answer.
`-M print` shows the atoms of each embedding. If a long query will not match, `-v -v -v` reports, every time
a match is attempted, the number of query atoms that were matched, which is usually a good
pointer to where the query fails. Note that query atoms are matched from left to right, following
branches depth first.

## Preparing molecules before searching

Several options change the molecule before it is searched

| Option | Effect |
| --- | --- |
| `-l` | reduce to the largest fragment first. A query for a counter ion will then not match |
| `-c` | remove chirality first. A query for `[C@H]` will then not match |
| `-g <name>` | apply a chemical standardisation, such as `-g nitro` to convert charge separated nitro groups to `N(=O)=O`. `-g help` lists them and `-g all` applies all |
| `-M omlf` | only keep matches that are in the largest fragment |
| `-X <element>` | remove atoms of that element |
| `-T E1=E2` | transform one element to another, see `-T help` |
| `-A <qualifier>` | aromaticity definition, see `-A help` and [aromaticity](../Molecule_Lib/aromaticity.md) |
| `-M imp2exp` | make implicit hydrogens explicit |

An example of `-g`. For the molecules `C[N+](=O)[O-]` and `CN(=O)=O`, the query `N(=O)=O`
matches only the second, and with `-g nitro` it matches both.

## Behaviour worth knowing

These are things that are easy to trip over. They describe the code as it is at the time
of writing.

* **Exit status does not signal a match.** The exit status is 0 whenever the run completes, whether or not any
  molecule matched. The non zero values are for errors, such as an invalid SMARTS or `-g`
  qualifier (61), a query file that cannot be read (6), or no input file (8).
* **Matches are counted as embeddings.** See [What counts as a match](#what-counts-as-a-match).
* **A single bond matches inside aromatic rings.** See
  [Aromatic bonds and Kekule forms](#aromatic-bonds-and-kekule-forms).
* **`-B` depends on the order of queries.** See [Several queries](#several-queries).
* **`-l` and `-c` change what can match.** A query for a counter ion or for a stereocentre
  will not match after the molecule has been reduced.

## Option summary

The options as listed by `tsubstructure`. Those that are not described above are for
specialised use, and the sub help they point to is the reference.

| Option | Meaning |
| --- | --- |
| `-s <smarts>` | SMARTS query |
| `-q <query>` | query file, with the prefixes shown in [Specifying queries](#specifying-queries) |
| `-m <stem>`, `-n <stem>` | write matches, non matches. `-m help` for qualifiers |
| `-o <type>` | output type, `smi`, `sdf` and others |
| `-i <type>` | input type and qualifiers. `-i help` |
| `-j <...>` | label matched atoms. `-j help` |
| `-J <tag>`, `-y <nbits>` | write fingerprints, and their size |
| `-a`, `-Y <stem>` | per query table, and its column names |
| `-P <type>` | atom typing used to determine changing atoms and match conditions (default `UST:AZUCORS`). `-P help` |
| `-f`, `-u`, `-k`, `-r` | what counts as an embedding |
| `-c`, `-l` | remove chirality, reduce to the largest fragment |
| `-M <...>` | miscellaneous query conditions. `-M help` |
| `-R <file>` | report file |
| `-b`, `-B` | stop testing queries early, see [Several queries](#several-queries) |
| `-G <file>` | write the matched atoms |
| `-A <qualifier>` | aromaticity. `-A help` |
| `-g <qualifier>` | chemical standardisations. `-g help` |
| `-X <element>` | remove atoms of this element before searching |
| `-T E1=E2` | element transformations. `-T help` |
| `-E <symbol>` | create an element with the given symbol |
| `-v` | verbose output |

Commonly used `-M` qualifiers, others are in `-M help`

| Qualifier | Meaning |
| --- | --- |
| `mmaq` | must match all queries |
| `kekule`, `nokekule`, `fkekule` | aromatic bond matching, see above |
| `minat=N`, `maxat=N` | only search molecules with at least or at most N atoms |
| `maxe=N` | at most N embeddings for each query |
| `stopm=N` | stop after N molecules have matched |
| `print` | print every embedding |
| `omlf` | only matches in the largest fragment |
| `imp2exp` | make implicit hydrogens explicit |
| `app=<text>` | append text to the names of the molecules that match |
| `rx=<regex>` | only process molecules whose name matches the regular expression |
| `report=N` | report progress every N molecules |
| `time` | report timing |
| `owdmm`, `anmatch` | adjust the `-a` output, see [Per query results](#per-query-results) |

## Searching large files in parallel

A search can be fast. A simple query may take about two seconds on 100k molecules, including reading
them. Larger files can be split and searched in parallel with `dopattern.sh`, or
tsubstructure_parallel.sh. See the [common tasks](../CommonTasks.md) page.
