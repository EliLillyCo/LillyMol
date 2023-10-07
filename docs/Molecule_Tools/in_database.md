# in_database

## Problem
Companion programmes `buildsmidb` and `in_database` build and use
key/value databases to enable structure matching by unique smiles.

Given a new set of molecules, one can then ask the question 'which of this
new set of molecules is identical to something in the database`. Frequently
the database might be a corporate collection, and the query molecules
a set of external molecules.

## Databases
Initially databases were built with gnu `gdbm`, and later migrated to
Berkeley DB databases. Given potential license difficulties with
Berkeley DB, there are efforts to identify a new key/value store
database that can serve as the underlying database. Leveldb has
been tried, but found to be wanting - interrupted operations seem
to leave the database in an unusable state, and concurrent reads
seem to be unsupported. Currently I am exploring lmdb, but may try
open source Scylla. For now Berkeley DB is the main target.

In all current use cases, the databases are built and then used in
read only mode. Something like the corporate database gets rebuilt once
per day, something like CHembl gets updated only when new Chembl releases
happen. Therefore there are no considerations of concurrent read/write
access.

# Uniqueness
Different measures of identity will be needed at different times. The
following attributes can be optionally set or unset when making structure
comparisons.

* presence of counterion(s)
* presence of chirality, @ or / \\
* presence of isotopes
* comparison as molecular graph

# Usage
The following options are recognised
```
Checks for molecules in a Berkeley database by looking up the unique smiles
Usage: in_database -d <dbname> <options> <input_file>...
  -d <dbname>    specify Berkeley database
  -j             discern lookup transformations from stored info (recommended)
  -H ...         specify storage conditions, enter '-H help' for info
  -b             look for both chiral and non-chiral smiles
  -l             strip to largest fragment
  -F <file>      stream for molecules found in the database
  -p             append database identifier to name
  -e <sep>       search all databases, insert <sep> string between database finds
  -U <file>      stream for molecules NOT found in the database
  -M <file>      stream for molecules with database marker
  -D <file>      write tabular output of search results
  -m <string>    specify database marker (default 'IN DATABASE')
  -W <fname>     write per atom count lookup rates to <fname>
  -n             ok to have no output - informational lookup only
  -Y             echo the DB key (useful for debugging)
  -i <type>      specify input file type
  -o <type>      specify output file type(s)
  -K             skip smiles errors
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -T ...          standard element transformation options
  -A <qualifier> Aromaticity, enter "-A help" for options
  -v             verbose output
```

## -d \<dbname\>
Specify one or more databases that have each been built using `buildsmidb`.
If there have been incompatible storage conditions specified during those
builds, an error will happen. By default, only the first database match
is returned, but see the `-e` option below.

## -j
when a database is built, the storage conditions are written to the
database. The most robust usage is to use those same storage conditions
during lookups. For example if a database was built, excluding chirality,
then chirality should be excluded during lookups.

## -H <conditions>
Manually specify storage conditions on the command line. This is risky
and mostly used for testing and debugging.

## -b
Look up both the chiral and non-chiral forms of the input molecules.
If a chiral match is found, that is returned. If that fails, chirality
is removed from the query molecule and the lookup repeated.

## -l
Strip to the largest fragment before doing a lookup.

## -F \<fname\>
Write the molecules found in the database to \<fname\>.

## -p
When writing the `-F` file, append the database identifier to the
molecule name. 

## -e \<sep\>
When multiple databases are present, search them all, and insert \<sep\>
between the database entries that match. Given input like
```
C carbon
```
the resulting output might look like
```
C carbon:methane:CH4
```
where two database matches have been found, one that called it `methane`
and the other called it `CH4`.

## -U \<fname\>
Write the molecules not found in the database to \<fname\>.

## -M \<fname\>
Do not use, the output seems redundant.

## -W \<fname\>
Sometimes it can be interesting to see how retrieval rates vary with
atom count. Accumulate that data and write to \<fname\>.

## -n
Normally `in_database` must be given an output file, either `-F` or `-U`.
But sometimes as a diagnostic, it can be useful to run with no output.

## -Y
Before each lookup, echo the key used for the lookup. This is useful
if transformations are being done on the smiles. Mostly for debugging.

## -T \<etrans\>
Apply element transformations prior to database lookup.

The `-O` option is an idea that really did not work. The idea was that one could
build an associated database of molecular formula. A molecular formula is much
cheaper to compute than a unique smiles, so the idea was to first lookup the
molecular formula, and if that was not found, then no need to do the expensive
unique smiles calculation. In retrospective, I don't think the molecular formula
is specific enough, and usually this makes thing worse.

## Building the Database
companion tool `buildsmidb_bdb` is used to build the BerkeleyDB databases
consumed by `in_database_bdb`. Many of the options are the same.
```
Usage: buildsmidb_bdb <options> <input_file>...
  -d <dbname>    specify BerkeleyDb database
  -S <fname>     write to a hash and then write hash
  -O <dbname>    optional database for molecular formulae
  -H ...         specify storage conditions, enter '-H help' for info
  -b             store chiral and non-chiral smiles
  -l             strip to largest fragment
  -i <type>      specify input file type
  -n <separator> separator for when storing duplicate entries
  -p             don't store duplicates
  -D <file>      write duplicates to <file>
  -U <file>      write non-duplicates to <file>
  -o <type>      specify output type(s) for -D and -U options
  -r <number>    report progress every <number> molecules
  -C <fname>     write changed keys to <fname>
  -x <natoms>    discard molecules containing more than <natoms> atoms
  -y             suppress messages about bad valences
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -T ...         standard element transformation options
  -A <qualifier> Aromaticity, enter "-A help" for options
  -v             verbose output
```

While the most common usage is to start with an empty database, this tool
creates an empty database if the `-d` file does not exist, it is also possible
to load into an existing database.

In both cases, decisions about how to handle duplicate structures must be made.
The default mode is to append the new name to the stored name, so the database
becomes a mapping from a unique smiles form, and one or more identifiers.

So, during a lookup, there may be a list of molecules returned as equivalents.

A minimal build will include the name of the database, and optionally some
combination of the structure modifiers. All other settings are optional. The
main thing is that there needs to be consistency between the options used during
building and lookup.

Molecular transformations are stored in the database (see below) and perhaps some
of this can be automated, but I think it will always involve some level of complexity.

### Performance
In the Lilly environment, we frequently use a globally mounted, NFS file system. While
this works well for many tasks, database building is generally better done in `/dev/shm` if
you have sufficient RAM available, or on a locally mounted scratch file system. For small
databases, this may not matter much, but for larger systems it does.

When loading large databases, performance on a the NFS system can become abysmal. The load
will start off seemingly OK, but as data is added, the process becomes slower and slower.

The tool `iwbdb_cat` can be used to concatenate BerkeleyDb databases. It can also operate
in a mode where duplicate values can be concatenated with existing values. The only
complication is that `buildsmidb_bdb` stores a special key, `_STORE_INFO` that contains
the molecular transformations required for matching. If multiple databases are
concatenated, that key will become unusable. One possibility is to allow the
concatenation to happen, and then manually fix the bad `_STORE_INFO` key at the end.

During lookups, `in_database_bdb` is fast. On a shared system, the more people use it
the better the file system caching becomes, and the faster results are returned. Unlike
what might happen with traditional database systems, where a server process can become
a bottleneck. One of the tests we frequently perform, the validate the integrity of
the process, is to lookup something like Chembl in a BerkeleyDb database of Chembl. That
usually takes under 5 minutes on very old hardware, and with an NFS mounted database.

Again, if performance is an issue, a locally mounted disk is best - or even /dev/shm
if you have the RAM to do it. We have often temporarily copied a BerkeleyDb database
to /dev/shm in order to do a specific large scale lookup.
