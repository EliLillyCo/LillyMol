# Ring Rarity
`smi2rings_bdb` is a tool for identifying new, rare, or well
precedented rings and ring systems.

Using tools like AI generative models to make molecules can lead to
unrealistic molecules. One of the most effective ways of identifying
unrealistic molecules is to determine whether they contain unprecedented,
or minimally precedented rings, or ring systems.

Note that this tool was developed several decades ago, for the inverse
problem: when purchasing molecules from external sources, molecules that
provided unique rings and ring systems were considered desirable. Today,
computer tools that generate unique rings are undesirable.

Note too that most existing collections will have considerable numbers of what can
best be called data entry problems, and these will likely include
ridiculous looking rings. See below about prevalence.

This tool works on ring systems only, so in this document, all references
to 'ring' should also include 'ring system'.

## HOWTO
In order to assess precedent for a new molecule, database(s) of already
precedented rings must
first be built. Prefer to use collections of molecules where there is
high confidence that the molecules have been made. Profiling collections
of virtual molecules raises the possibility of introducing rings
that someone thought could be made, but in fact it was impossible.
A corporate collection would be a good candidate, as would Chembl molecules
where there is associated assay data. Etc...

There is of course no 'right' answer to what collections to consider. Perhaps
in a generative study, if people have speculated that a ring is possible,
maybe that is good enough - since this test will identify rings that
nobody has ever contemplated. But why has that ring not appeared in 
collections of actual molecules?
This is complex and you will need to make decisions.

### Granularity
What is a ring? This tool provides three descriptions of a ring - or
ring system.

1. The ring with zero information about substitution.
2. The ring with attachment points marked
3. The ring with attachment points marked by what kind of atom is attached.

There is no 'right' answer to which is most useful, it will depend
on the situation. All have proven to be useful at various times. You will
need to decide what resolution makes sense for your situation.

Of course if a new molecule contains a ring, where examples of that ring
with exactly the same kinds of substituents have been observed, there is high
probability this molecule could be made. If the ring has been observed,
but never with this substitution pattern, is that OK or not? Again,
there are no 'right' answers.

### Build
When building a rings database, extract all representations, which one
to use can be decided at lookup time.
```
smi2rings_bdb -d STORE -d /dev/shm/collection.rings.bdb -g all -v -n \
    -Y hash -j ring -j iso -j env=UST:achry -j double -j spiro \
    -z 8 -N add collection.smi
```

### -d STORE -d /dev/shm/collection.rings.bdb
Specify the database to be built. Note that large Berkeley databases can
suffer performance problems if built on nfs mounted file systems.

### -g all -v -n
Chemical standardisation, verbose, and suppress normal standard output -
which is a list of the rings produced.

### -j options

-j double: Append exocyclic doubly bonded atoms to the ring. This will be very
important for aromaticity. This should be used all the time, perhaps it should
not be an option.

-j spiro: Ring systems span spiro fusions - by default in LillyMol ring systems
do not span spiro fusions. In this situation you probably do want this.

-j iso: Generate the version with an isotope (1) at the attachment point(s). During
lookup this means that a molecule will only match if all the ring atoms are the
same, but the substitution pattern must also match.

-j ring: store/retrieve the unsubstitued ring. No information about substitution.
During lookups, a ring will match regardless of the substitution pattern.

-j env=UST:achry: Generate a version where the atom type of the attached
atom is placed as an isotopic label in the ring. This will likely lead to
ridiculous looking isotopic labels, `N1C[12015CH2][4525CH2]CC1` for example.

The atom type to use, 'achry' in this case, can be any atom typing supported
in LillyMol. There is no 'right' or 'best' atom type to use. This one includes
the aromaticity, number of connections, number of hydrogens, ring and compressed
atomic number.

-z 8: do not store any ring with more than 8 atoms. Note that this means that
during lookup if a ring contains more than 8 atoms, it will not be found.

-N add: In the case of [n+] atoms, add the connected atom. You probably want this.

-Y hash
By default, smi2rings_bdb writes each ring to the database as it is found.
If that ring is found again, the number of examples already encountered
must be fetched from the database, the count incremented and the proto is
written back to the database. That can be quite inefficient. The `-Y hash`
combination holds all ring information in memory and only writes the
database once all molecules have been read.

For example processing a collection of about 2.7M molecules 6:31
minutes with the default invocation but only 5:09 using the hash.

-Y scaffold
The ring systems of a molecule are just a well defined subset of atoms.
The scaffold is another such well defined set, so instead of storing
ring systems, the scaffold can be stored and used for lookups.
Note that if you build a database with `-Y scaffold` and then want
to do scaffold lookups, you must also specify that during lookups
```
smi2rings_bdb -d /path/to/db.bdb -Y scaffold -d LOOKUP ... unknown.smi
```
It is an interesting thought experiment to consider what would happen if we
stored the scaffold complement, the molecular 'spinach' in the database
instead. Given an active molecule we could then find other molecules
that presented the same functional groups, but with different
cores. Kind of like ring replacement, but for scaffolds.

Other options to consider:

-c remove chirality.

This might depend on what you believe to be the reliability of chirality
information available, and whether or not a generative model is making
chiral molecules. I should probably implement the ability to store
both chiral and non-chiral forms...

#### Building databases

Since the identity of the first molecule exemplifying a given ring is 
stored in the database, it is usually desirable to sort the collection
by heavy atom count and use the sorted file. Use `msort`.

For each collection, build a rings database. Note that the data in the
database is a textproto representation, the unique smiles of the ring being the key.
Depending on the -j options used during the build, this may be an un-annotated
ring (-j ring), a ring with isotope 1 at the attachment points (-j iso), or isotopes
indicating the atom types of what is attached (-j env=...).
```
C1=C2C3CCCCC3CCC2CCC1 n: 25 ex: "CHEMBL116328" 
C1[1CH2]CC[1CH2]1 n: 3 ex: "CHEMBL168248"
[3518cH]1ccccc1 n: 77264 ex: "CHEMBL501101"
```
with the seemingly strange isotope, 3518, being an atom type.

Use `iwbdb_list` to
examine the contents of a Berkeley DB database. The resulting output is
a tabular text file and can be processed with many utilities. And since
the first column is a smiles, structure aware tools can be used. Beware however
that there will be instances of aromatic forms that cannot be decoded.

When operating in
lookup mode, `smi2rings_bdb` can open multiple databases, just specify
multiple -d options. But only one database during building.

### Concatenation
As we see above, the count is the second token in the database value. So
if you wish to concatenate multiple databases it can be done with
```
iwbdb_cat -d result.bdb -v -R 2 split1.bdb split2.bdb ...
```
Or if you have multiple collections, use this to make an aggregate
database, to avoid having to open multiple databases. Note that
using the column ordering is not robust, better would be to
write a tool that parsed the proto.

During lookup, it can be handy to know how the processing is progressing, use
`-Y rpt=<nn>` to have smi2rings_bdb report every <nn> molecules processed.

## Lookup
Once the known collections have been profiled, new molecules can be looked up
```
smi2rings_bdb -d LOOKUP -d collection1.bdb -d collection2.bdb... -v unknown.smi
```
specifying via the -j option(s) which ring forms to lookup

### Output
Apologies in advance, the default output is kind of horrible. We have just added
a textproto output format which should make parsing easier.

Current default output might look like
```
O1CC(N=C1C1=NC(=CC=C1)C1=NC(C(C)C)CO1)C(C)C PBCHM500120 URS:3 CHEMBL3652777 3, 6
```
The first tokens are the smiles and ID of the input molecule.

Following is a URS:<n> token that us the number of instances of the rarest
ring in PBCHM500120, in this case in that rarest ring was found
in the Chembl rings database, id CHEMBL3652777.

If the ring is found, the name of the exemplar structure in the database will
be reported.

Following that will be a variable number of tokens, according to the number
of ring systems in the input molecule, and the 'URS' score for each one.

Note that in practice most databases are built with a limit on the ring size
stored. This can create confusion if a query molecule with a ring larger than what
is stored in the database is used. The output will show 'URS:0' which means
the molecule contains a ring not in the database. But that molecule may have
been in the file from which the database was built. It may be convenient to
add `-Y applarge` which appends to each molecule with URS:0 the size of the
largest ring in that molecule.

#### -Y proto
This enables output of Smi2Rings::Results protos.
```
smiles: "N1(C)N(C)CC1N" id: "PBCHM134269649" result { urs: 64 ex: "PBCHM323358" }
```
In this case there are 64 examples of this ring in the database, with
PBCHM323358 being an exemplar of that ring. If the '-a' option is not
given, there will be a result for each ring system. These will be
sorted by rarity.

#### -Y smiproto
This writes a regular smiles file with the Smi2Rings::Results proto appended
to the molecule name
```
N1NC(N)CC1C PBCHM134081637 result { urs: 28305 ex: "PBCHM4867587" }
```
In this case the ring in that molecule is exemplified by PBCHM4867587
and there are 28k instances in the database.

### Prevalence
You can use the tool as a filter with the '-f' option. This allows specifying
a minimum 'URS' value for output. Note again that most collections contain
likely data entry errors. In practice, we have found a value of 5 or 10 for '-f' to
be a reasonable compromise. Adjust to taste.
