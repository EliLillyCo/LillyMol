# Firstatom

# Purpose
Firstatom write smiles where an atom, specified by a substructure
query, is the first atom in the smiles. The most common use case
for this is where there might be disparate functional groups across
a series of molecules, but a particular one must be processed.
`firstatom` can be a preprocessing step, possibly quite complex,
that places the right atoms first in the smiles.

## HOWTO
The usage message is
```
Molecule_Tools/firstatom.cc compiled redacted redacted
Creates smiles starting with an atom hit by a query
  -q <query>     specify query or queries
  -s <smarts>    specify one or more smarts queries
  -m <number>    query atom number to use (default 0)
  -z i           ignore molecules not hitting any queries
  -z w           write molecules that don't hit any queries
  -z f           take the first of multiple hits
  -z m           skip queries with multiple hits
  -z each        when a query matches multiple times, do each match
  -k             do NOT perceive symmetry related matches
  -u             unique embeddings only
  -P <prefix>    prepend <  prefix> to all smiles (useful for adding special atoms)
  -b             write molecular subset (matched atoms only)
  -i <type>      specify input file type. Enter '-i help' for details
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -A <qualifier> Aromaticity, enter "-A help" for options
  -v             verbose output
```

Usage is pretty straightforward, there is not much to explain.

## Options
* -q \<query\>   a query specification.
* -s \<smarts\>  the query(s) entered as smarts
* -z i           ignore molecules not matching any of the queries.
* -z w           write non-matching molecules to the output
* -z f           when there are multiple query matches, arbitrarily write the first match.
* -z m           skip molecules that have multiple matches.
* -z each        when multiple matches present, write each match as a separate molecule.
* -k             do not perceive symmetry related matches
* -u             unique embeddings only.
* -P \<prefix\>  prepend \<prefix\> to all smiles - this is a text addition.
* -b             write the matched atoms subset, prefer to use molecule_subset.
