# Welcome to the Eli Lilly LillyMol implementation.

## Background
LillyMol is a set of Linux executables for Cheminformatics. These tools are built
on a high performance C++ library for Cheminformatics.

LillyMol does only a subset of Cheminformatics tasks, but tries to do those tasks
efficiently and correctly.

LillyMol has some novel approaches to substructure searching, reaction enumeration and
chemical similarity. These have been developed over many years, driven by the needs
of Computational and Medicinal Chemists at Lilly and elsewhere.

Recent work has focussed on making *de-novo* molecule construction and there are
several tools desiged to either support or complement A/I driven molecule
generation.

LillyMol is fast and scalable, with modest memory requirements.

This release includes a number of C++ unit tests. All
tests can be run with address sanitizer, with no problems reported.

The file [Molecule_Tools/introduction.cc](src/Molecule_Tools/introduction.cc) provides
an introduction to LillyMol for anyone wishing to develop with C++.

This release includes first steps towards more extensive documentation of LillyMol,
see the [docs](/docs) directory. More work is needed on this front. Most parts
of LillyMol have been feature stable for a long time.

## Python
This release include a python interface to LillyMol via pybind11. This first release
includes most Molecule related functionality, substructure searching and reaction
enumeration. In the pybind directory there are some *_test.py files that exemplify
much of the current functionality. Documentation is in [docs](/docs/python/LillyMolPython.md).
This should be sufficient support for a great many tasks involving querying
or manipulation of molecules at the connection table level.

The current roadmap for the python interface primarily involves two directions

* Enabling gfp fingerprints for similarity calculations.
* Making existing LillyMol applications available.

There is already a Julia interface to an earlier version of LillyMol, and this
release will soon be adapted to support Julia.

See [BUILD](Build.md) for how to build LillyMol.

## Structure Files
All LillyMol tools make certain assumptions about their input file(s).

Smiles files are assumed to look like
```
C methane
CC ethane
```
where the smiles and the name are separated by whitespace - space or tab.
Generally LillyMol supports spaces as the preferred delimiter.
Multiple spaces/tabs between smiles and id are ignored.
Smiles files to not have header records.

All tokens after the smiles are part of the molecule name.
```
C methane consists of           a Carbon atom
CC methane has two.
```
is a Methane molecule whose name consists of 6 words, and an Ethane whose
name consists of 3 words.
Note that spaces in the name are NOT compressed, so the name of Methane
above is quite long - containing a lot of whitespace.
Beware however that if you
consume these molecules and write out a result that includes the name,
this may be quite problematic - some tools will write only the first token
of the name, some may change the spaces to underscores, some may die.
To be safe, work with single token molecule names.

### SDF Files
See [IO](/docs/Molecule_Lib/io.md) for more information about reading .sdf files
with LillyMol.

### Tabular Files
LillyMol descriptor generation tools generate tabular files that typically look like
```
Id D1 D2 D3
id1 1 2 3
id2 4 5 6
id3 7 8 9
```
which are space delimited and with a header record. Use tcount to check whether
or not a file is tabular.

## Functionality

When encountering a new set of molecules, the first task is often to transfer from
a format like sdf. For that fileconv can be used.
```
fileconv -i SDFID:cat_num -f lod -E autocreate -e -O def -I 0 -V -g all -c 6 -C 50 -v external.sdf
```
this will create 'external.smi' with the results of the conversion. We specified
that the identifier we care about is the sdf tag 'cat_num'. 

Using '-E autocreate -e' says if atoms like `[R]` or `[Pol]' or `[X]` or .... are
encountered, do **not** die, but consume them as valid molecules (-E autocreate) but then
do not write them (-e).

We do not want isotopes, so `-I 0` does that. We could of course convert them with `-I change`.

The `-f lod` directive says reduce these molecules to the likely largest fragment of interest.

We alomost certainly do not want molecules with valence errors, `-V` discards those.

We apply [chemical_standardisation](../docs/Molecule_Lib/chemical_standardisation.md) via the
`-g all` option.

Discard molecules with fewer than 6 or more than 50 atoms. Note that `-g all` will remove
any explicit Hydrogen atoms, so the atom count cutoffs apply to heavy atom counts.

If you have a large number of molecules, it can be helpful to use `fileconv_parallel.sh`
but in this case you must specify the name of the output file
```
fileconv_parallel.sh <same as above> -thr 16 -S external external.sdf
```
will run 16 individual fileconv tasks on your computer - make sure you actually
have 16 cores!

## Molecular Properties
In order to understand what is in a collection of new molecules, generating
molecular property profiles can be very informative 
[molecular_property_profile](/contrib/Molecular_Property_Profile/README.md).
This shows aggregate collection wide properties such as distributions of
several interpretable molecular properties.

## Overlap with Corporate
Usually one of the most interesting things about a new set of molecules is
how does it overlap with the existing corporate, or public, collection(s).
Assuming you have used [buildsmidb_bdb](/docs/Molecule_Tools/in_database.md)
to build a database of those existing collections, use in_database_bdb
to look up the new molecules in those databases.
```
in_database_bdb -d /path/to/corporate/all.bdb -d /path/to/chembl.bdb \
        -F in_database -p -U notindatabase ... external.smi
```
where you will need to use the same options during lookup as were used
when building the database - dropping chirality, fragment selection, tautomer
matches, etc...

## Similarity with Corporate
For each new molecule, how close is it to something you already own?
Generate fingerprints
```
gfp_make.sh external.smi > external.gfp
```
The tool `gfp_lnearneighbours` is used to compare two sets of fingerprints.
Indeed it can be used to compre two large collections, but it will take
forever. In this case split up the external file into chunks of (say)
5000 fingerprints
```
iwsplit -n 5000 -tdt -suffix gfp external.gfp
```
which will generate a bunch of `iwsplit*.gfp` sharded files containing
everything from 'external.gfp`. At this stage some kind of parallel
processing is needed. Inside Lilly it is common to use the SGE
cluster for this task. That might look like (assuming 237 splits formed)
```
dopattern.sh -cluster -o 237 'gfp_lnearneighbours -n 1 -p iwsplit%.gfp \
                /path/to/corporate/all.gfp.gz > iwsplit%.nn'
```
which will generate a bunch of iwsplit*.nn files, each one containing
the nearest neighbours from a portion of external.gfp. These files
can be combined
```
nplotnn -H histogram.txt -v iwsplit*.nn > all.nn
```
That file can be sorted
```
tdt_sort -T DIST all.nn > all.sorted.nn
```
so the fingerprints with closest neighbours are at the top of
the file.

The file creaed with the -H option contains a histogram of nearest
neighbour distances which can be plotted to get an idea of how
similar are the new molecules compared to existing collections.

[!NOTE]
It is quite possible for the external collection to be quite different
from the existing collection, but have very low internal diversity. This
measure is only about how different from the existing collection.

## Internal Diversity
The internal diversity of a set of molecules is hard to assess.

The first thing to check is for duplicate structures.
```
unique_molecules -S unique -c -z -l -v -t I=Cl -t Br=Cl external.smi
```
We are not interested in chiral variants, chirality and cis-trans bonding is discarded, -c -z.
Counterions are discarded (if not already happened) with the -l option. We don't
really consider the Iodo and Bromo variants of a Chloro group to be "differen
enough" so we translate all the heavy halogens to Chlorine. Again
there is unique_molecules_parallel.sh for larger sets of molecules.

But exact matches are only part of the uniqueness story. Similarity
based methods can be used.

If the set of molecules is not too large, say 2M or so, then use gfp_spread.

First generate fingerprints
```
gfp_make.sh external.smi > external.gfp
```
then run spread with as many cores as you have available, but no more than 16.
```
gfp_spread_standard -h 8 external.gfp > external.spr
```
and then plot the trajectory of distances [gfp_spread](/docs/GFP/gfp_spread.md).

### Unique Ring Systems
If you have built a database of the rings in your corporate collection, you
might be curious to know if there are novel ring systems in this new set
of molecules
```
smi2rings_bdb ... -d /path/to/corporate/rings.bdb -v external.smi > rings.smi
```
see [smi2rings_bdb](/docs/Molecule_Tools/smi2rings.md) since this will
require concordance with how the database is built and how it is queried.

You might be curious how many different ring systems are present
```
molecular_abstraction -a 'scaffold(WRITE)' external.smi \
        | unique_rows -c 1 -v - > unique_scaffolds.smi
```
unique_rows looks at the contents of column 1, the smiles, and discards 
duplicates. You can then come up with a term we have called the "scaffold
density" which is the total number of molecules, divided by the number of
unique scaffolds, or for each scaffold, on average, how many molecules
were made? The -O option of unique_rows will yield a table of scaffold
and number of occurrences.

## Property Filtering
We might have ideas about desirable ranges of molecular properties and wish
to discard molecules falling outside those ranges. Use molecule_filter,
or molecule_filter_parallel. Create a textproto describing the properties
we want, with min/and or maximum values.
```
min_natoms: 15
max_natoms: 30
min_heteroatom_count: 3
min_heteroatom_fraction: 0.2
max_heteroatom_fraction: 0.8
min_nrings: 1
max_nrings: 4
min_aromatic_ring_count: 1
max_aromatic_ring_count: 3
min_tpsa: 70
min_alogp: 0.5
largest_ring_size: 7
max_ring_system_size: 3
exclude_non_organic: true
exclude_isotopes: true
max_chiral: 2
max_halogen_count: 5
max_distance: 25
max_hbd: 10
```
Then invoke molecule_filter_parallel with this input file
```
molecule_filter_parallel.sh -F /path/to/mf.textproto -thr 16 external.smi > filtered.smi
```
takes 6.2 seconds to process 2.8M molecules. This of course only works if you
actually have 16 cores. If not, split the input file into chunks with iwsplit
and use dopattern to submit to the cluster.
```
iwsplit -suffix smi -n 5000000 external.smi
dopattern.sh -o 200 -cluster 'molecule_filter -F /path/to/mf.textproto iwsplit%.smi > filtered%.smi'
```
assuming 'external.smi' was split into 200 chunks each with 5M smiles. Seems like a lot...

While molecule_filter has only a relatively small number of properties upon
which filtering can be done, any of the properties computed by [iwdescr.sh](/docs/Molecule_Tools/iwdescr.md)
can be used as the basis for a filter. For example to filter on the size of the
largest electron rich section of a molecule
```
iwdescr.sh -F 'w_lercsct.gt.6' file.smi > passed.smi
```
And of course there can be more complex queries
```
iwdescr.sh -F 'w_htroaf.gt.0.1' -F 'w_aromdens.lt.0.5' file.smi > passed.smi
```
which filters to molecules where the heteroatom fraction is greater than 0.1
and the aromatic density is less than 0.5.

In both these cases these descriptors are in the default set, so adding `-O none`
will substantially speed the computation.


## Substructures Searching
When closely examining larger collections it is usually convenient to use
a sorted list of molecules
```
msort_parallsl external.smi > external.sorted.smi
``
2.3M Chembl molecules can be processed in 3.8 seconds.

There may have been a claim that this set of molecules contains aromatic ring systems
with both aniline and Halogen substituents, and no hydroxy groups. Build a query file
```
name: "aniline and Halogen"
query {
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      environment: "1a-[NH2]&&1a-[Br,I]&&0c-[OH]"
    }
  }
}
```
Note that the query specifies that on the ring system there should be exactly
one aniline group and exactly one of the Halogens, and zero occurrences of a Phenol.
We also restricted attention
to ring systems of size 2 which have two aromatic rings, so no aliphatic rings.

Then run the query
```
tsubstructure -q PROTO:/path/to/aniline_halogen.qry -m aniline_halogen external.sorted.smi
```

or use `tsubstructure_parallel.sh` if the file is large. This might match molecules like
![CHEMBL27011](Images/CHEMBL27011.png)
![CHEMBL81588](Images/CHEMBL81588.png)

And of course if a query can be specified as smarts, that can be used with the `-s` option.

We are looking for molecules that contain a likely primary amine '[NH2]-[CX4T1]' more
than five bonds away '...{>5}' from an acid 'C(=O)-[OH]'

```
tsubstructure -j 1 -s '[NH2]-[CX4T1]...{>5}C(=O)-[OH]' -m matches external.smi
```
which might match molecules line
![CHEMBL2286788](CHEMBL2286788.png)

Because we use the -j option to tsubstructure, the matched atoms are given isotopic
labels. Somewhat confusingly matched atom 0, the Nitrogen atom, is assigned isotope
1 - because isotope 0 is not an isotopic atom. Generally LillyMol uses zero based
indexing for matched atoms.


Or we want to mark a reactive Nitrogen atom in order or likely reactivity
```
tsubstructure -s '[ND1H2]-[CX4]||[ND2H]([CX4])[CX4]||[ND1H2]-a||[ND2H]([CX4])-a||[ND2H](a)-a' \
        -m amide -j 1 -v reagents.smi
```
The OR type queries in LillyMol, || are evaluated left to right. When one matches
matching stops, and the matched atoms will be returned and labelled with an
isotope. This might match molecules like
![CHEMBL4101550](Images/CHEMBL4101550.png
This match is interesting because it shows how the precedence matching has worked.
The primary amine is a match, and that was detected first
![CHEMBL70445](Images/CHEMBL70445.png).

## Reactions
When performing reactions, what is usually the hardest part is selecting the
reagents. `tsubstructure` allows for selecting molecules with just one functional
group. For example to identify molecules with just one acid
```
tsubstructure -s '1[OH]-C=O' -m acids -v all.smi
```
But beware, if you wanted molecules with precisely one secondary amind, this
will not work
```
tsubstructure -s '1[ND2H]([CX4])[CX4]' -m secondary_amines -v all.smi
```
Because the query is symmetric, this will always generate two matches. Either
change the numeric requirement to 2, or add the `-u` option to find only
unique matches.

And to avoid molecules with primary amines
```
tsubstructure -u -s '1[ND2H]([CX4])[CX4]&&0[ND1H2]-[CT1X4]' -m secondary_amines -v all.smi
```
or run sequential queries.
```
tsubstructure -u -s '1[ND2H]([CX4])[CX4]&&0[ND1H2]-[CT1X4]' -m - all.smi | \
     tsubstructure -s '[ND1H2]-[CT1X4]' -n secondary_amines -v -
```
The first invocation makes a positive match for a secondary amine and the
second only writes molecules that do **not** contain a primary amine. It is
upredictable which will be faster - on this system both took about 4.2 seconds
to process 200k random molecules. The shell is your friend!

As is often the case, there is often no 'right' way of doing a given task,
there may be several different approaches, each with advantages and
disadvantages.

Once the acids and amines are assembled, we might choose to perform
the reaction. That might look like
```
name: "acid amine"
scaffold {
  id: 0
  smarts: "[OH]-C=O"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "[1NH2D1]-[CC2T1]"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
Remember the amine had been isotopically labelled by tsubstructure, so
using it in the reaction becomes straightforward. I have seen a great many
truly incomprensible smarts in reactions written by very clever people trying to express
complex concepts. That can often be avoided via tricks like this.
```
tsubstructure -j 1 <very complex query> -m R1 start.smi
trxn -P acid_amine.rxn acid.smi R1.smi
```

Of course it can be done in one step if needed by using the `query_file` directive
rather than `smarts` in the reaction file. See [trxn](/docs/Molecule_Tools/trxn.md).
Using that construct, the exact same query that is used to filter the reagents
can be used for the reaction.

## Building a Model
We get activity data for a subset of our molecules in the file `train.activity`. That
file is space separated with a header
```
ID Activity
id1 3.14
id2 1.59
...
```
Fetch
those smiles from the larger smiles file
```
fetch_smiles_quick -v train.activity all.smi > train.smi
```
We now have train.smi and train.activity that have the same identifiers in them. Make
sure both are rectangular
```
tcount train.activity train.smi
```

Use [activity_consistency](/docs/Molecule_Tools/activity_consistency.md)
to explore the nature of duplicates, and possibly generate a more internally consistent
training set file.

Generate descriptors - restrict to some 2D descriptors
```
make_descriptors.sh -j 4 -abr -ap -hpo -w train.smi > train.dat
```
Remove constant and near constant columns
```
notenoughvariance -j -n 10000 -x 0.99 -v train.dat > train.nev.dat
```
Build an XGBoost model
```
xgbd_make.sh -mdir MODEL -A train.activity train.dat
```
Evaluate that model
``
xgbd_evaluate.sh -mdir MODEL -smi test.smi > test.pred
```
When given the -smi option, it assumes that the input smiles can be converted
to descriptors via make_descriptors.sh.

See how good the prediction is
```
iwstats -E test.experimental -p 2 test.pred
```

## Sphere Exclusion (leader).

Having scored a set of unknown molecules with our model, we wish to make 
a diverse selection from that set. Assign a relative desirability to each
molecule and add that as an extra column to the smiles file
```
smiles1 id1 score1
smiles2 id2 score2
.
```
There are an unlimited number of ways this could be done.

sort that file, most desirable at the top of the file, and generate
fingerprints on the sorted file.
```
sort -k 3,3g test.smi > test.sorted.smi
gfp_make.sh test.sorted.smi > test.sorted.gfp
```
That fingerprint file has the most desirable molecules at the top of the file.
gfp_leader implements a sphere exclusion selection where the first molecule
selected is the first molecule in the file. A fixed radius is specified,
and all unselected molecules that are closer than that distance to the first
selected molecule, are placed in that first cluster. Once the first cluster
is formed, the search returns to the top of the list and the first unselected
molecule becomes the next cluster centre (leader). Etc.

This way, we generate a set of molecules that are both desirable and diverse. We
can guarantee that the cluster centres are all more than the distance threshold
separated. And because the cluster centres have all come from the top of the
file, they are enriched in desirable molecules.

There is no 'right' radius to use. It depends on how much you trust the ranking
function, how many experimental slots you have, many things...

```
gfp_leader_standard -h 8 -t 0.20 test.sorted.gfp > test.sorted.ldr
nplotnn -n 0 -v -L def -v test.sorted.ldr > test.sorted.ldr.smi
```
If you only need to select a certain number of clusters (usually), add the
`-n` option to gfp_leader_standard. If you have a very large dataset
consider discarding all molecules that are below a certain level in
the sorted file. There is little chance that they could ever be
selected. Selecting 1000 molecules from a collection of 2.8M and
radius 0.20, with 8 threads, takes 25 seconds. With 16 threads that
drops to 19 seconds. While this time is OK, anything much larger
would benefit from discarding the least desirable candidates.

## Other Properties
# QED
```
QED.sh file.smi > file.qed
```
Alogp
```
alogp -M 6 file.smi > file.alogp
```
Rotatable Bonds
```
rotatable_bonds file.smi > file.rotb
```

Note that both alogp and rotatable bonds are computed as part of iwdescr,
although in the case of rotatable bonds, the computation is slightly different.

We notice that our new collection of molecules seems to have some very elongated
molecules. Extract some of those.
```
long_molecules -F long -d 3 -D 4 -m 20 all.smi
```
We are looking for molecules where the longest path is at least 20
bonds, the average distance from that longest path is 3 and the longest
distance from that path is 4 atoms. This matches molecules like
![CHEMBL13003](Images/CHEMBL13003.png)
We may wish to filter molecules like this - they might not fail other
filters, although a rotatable bond filter may get some of them.

## Geeks Only!
When working at low levels with molecules it is often convenient to
deal with atom numbers. The tool `numbered_smiles` consumes a connection
table and writes it with whatever is specified as the isotopic label.
So given the smiles
```
ClC(CCl)CO CHEMBL1538584
```
it will generate
![CHEMBL1538585](Images/CHEMBL1538585.png)
A number of other properties can be applied as isotopic labels.
