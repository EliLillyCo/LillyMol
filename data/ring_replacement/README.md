This directory contains a set of ring replacements extracted from a version
of Chembl.

The files are named according to the kind of ring or ring system. For example
'rings.5a.smi' contains a list of five membered aromatic rings. On the other
hand the file `rings.5A.smi` contains a list of five membered aliphatic rings
extracted from Chembl. This naming convention follows the aromatic/aliphatic
conventions used in smarts.

There are case insensitive file systems where file names that differ only in
case cannot coexist.

For this reason, the files in the distribution are "hidden" (start with .) and
need to be converted into the file names suitable for the current file
system.

There are two scripts

to_linux.sh
to_mac.sh

where the first one is for file systems where 5a5a can coexist with 5A5A
and the second is for case insensitive file systems. In that case, 5a5A will
show up as 5Ar5Al (Aromatic and Aliphatic) a much less pleasing naming scheme.
