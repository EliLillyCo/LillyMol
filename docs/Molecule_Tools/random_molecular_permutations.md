# Random Molecular Permutations

**Note documentation is work in progress**

Random molecular permutations is a chemistry aware tool that makes
random changes to input molecules. In terms of generative functions,
it is good at making local variations on starting molecules. 
The objective is to create new molecules that are plausible
variants of the starting molecules.

It contains a set of pre-defined molecular changes that are randomly
applied to the current molecule.

At each step, the tool chooses among a range of possible
transformations and attempts to apply that transformation to the
current candidate molecule.  The probability of selecting each kind of
transformation can be controlled via configuration options.

Note that by default, the tool does not have the ability to
favor a given outcome, but filters can be applied at each step
and undesirable molecules dropped. So while there will not have
been optimization of an objective, undesirable molecules
should be diminished. Fundamentally, the tool is dumb, but
hopefully useful.

## Transformations.
Currently there are about 25 different molecular transformations
defined. I am in the process of transitioning `random_molecular_permutations`
to proto configuration. The transformations will be described in the
proto, rather than here. The proto is
`Molecule_Tools/random_molecular_permutations.proto`.
The `TransformationProbabilities` message describes the transformations
available.

At each step, a transformation is randomly
selected. Then, based on its associated probability, it will
choose whether or not to change the molecule. In this way we
can make certain transformations rare - break a ring, while having
other transformations more common - add a carbon atom. Transformations
can be turned off by assigning them a zero probability.

Central to the functionality is the idea of external fragments. While
the tool can add individual C, N and O atoms, larger changes can
come from adding fragments from an externally defined fragment library.

## Libraries
As described in the `Libraries` message in the proto, several
fragment libraries are also used. These can be any external
file satisfying the requirements for a particular transformation.
Note that `random_molecular_permutations` has in-built transformations
that result in the addition of C, N and O atoms, so it will never be
necessary to include those as library fragments.

It is important to note however that the probability of any particular
fragment being added, can be quite low. First the add fragment from
library transformation will need to be selected. Then that transformation
will random choose a library member to use. Only if that library
member can be attached is a successful change observed. This could
be addressed by having multiple, independent fragment libraries,
but today, all libraries on the command line are concatenated into one.

## HOWTO
The tool is necessarily complex, due to the fact that there are
a large number of different kinds of transformations possible,
each of which has a probability.

There are two parameters that are fundamental to the operation of
`random_molecular_permutations`

1. The number of copies of the starting molecule, `-x` option
2. The number of recursive permutations of each copy to make `-p` option.

For each input molecule, sequentially make `-x` copies of
that molecule. Then for each copy make a single, randomly selected, change.
If the `-p` option is greater than 1, make random change
to that product molecule, and then if the `-p` option is greater
than 2 make...

Clearly as the `-p` option increases, molecules will diverge
from their starting point.

Not all attempts at making a permutation will be successful. The
molecule, either the original, or one of the permuted variants,
might not have atoms suitable for a specific change, so some
degree of inefficiency is inevitable.

## Filtering
Molecules generated are checked to see if they violate
certain rules that are set. These can be specified on
the command line

* `-c` lowest number of atoms in a molecule
* `-C` largest number of atoms in a molecule
* `-r` lowest ring size allowed
* `-R` largest ring size allowed
* `-y` lowest number of rings allowed
* `-Y` largest number of rings allowed

### Substructure Based
Substructure queries can be used as either positive or negative
filters.

If the `-q` option is present, only molecules that match **all** of
those queries will be written or passed on. I am thinking this
`all` match should be an or match...

If the `-Q` option is present, any molecule that matches one of
those queries will be immediately discarded.

Note that if you want to specify these queries as smarts, the syntax
`-q 'SMARTS:[CD3H](=O)-[#6]'` works.

Note too that *or* conditions can be specified with smarts like
```
-q 'SMARTS:O=[CD3]-N-a||COc'
```
which looks for an aromatic amide, or an aromatic methoxy. Which could
also have been done via a much more complex recursive smarts.

## Objectives
One of the main problems with this tool is that it does not
grow towards any objective, molecular variants are generated
and the results are not fed back into the generation process.
To some extent that is OK, because it is unclear that if
removing a CH3 group were good, then removing another one would
be similarly advantageous. But if removing a CH3 from one starting
molecule were advantageous, then it might be advantageous for another
starting molecule.

It is not all bad, if filters have been specified, molecules violating those
filters do not progress further. So, there is some mid-course direction,
but only via negative selection, with nothing to guide molecule
generation in a desired direction.

### Similarity Based
Use the -T option to specify a set of target molecules to aim for. This is
the only positive selection currently implemented. It is complex. The usage
message for this is
```
Options for controlling building towards a target population of molecules
 -T <fname>          file containing the desired molecules
 -T nmax=<n>         the number of proximal molecules to retain
 -T dmax=<d>         only retain proximal molecules that are closer than <d>
 -T fdbk=<n>         feed back newly identified proximal molecules <n> times
 -T tom=<d>          turn off any target molecule matched to within <d>
 -T write=<fname>    write all molecules added to the proximal collection to <fname>
```

A set of molecules is read from <fname> and their fingerprints are computed.
As new molecule are generated, they are compared with this target set of molecules.
Molecules that are further than <dmax> are ignored, but still continue
being processed.

The `tom` directive means that once a target molecule has been matched,
within `dmax` from a newly generated molecule, it is deactivated. This might
help avoid saturating the proximal list with molecules close to
a single entry in the target set - see `nmax`.

If the `nmax` directive is used, a list of newly generated molecules,
that are closer than `dmax` from the target molecules is maintained. 
This list always contains the `nmax` molecules closest to a target
molecule found so far.

Likely you will want to write these molecules via the `write=` directive.

If the `fdbk` directive is given, these proximal molecules are started
as their own seeds, which should help explore molecules close to the
target set.

A typical invocation might be
```
-T desirable.smi -T dmax=0.15 -T nmax=200 -T fdbk=20 -T write=newly_found_desirable
```

# Future
Clearly the generation and filtering process could be enhanced with tools such as

* ring rarity
* synthetic precedent
* `w` descriptors
* lookup in database.

all of which are `C++` and could be linked.

I always wanted to hook up the tool to a model, and use the model feedback to
decide whether a proposed change was improving things or not. Allow for a
few mis-steps, like in a simulated annealing type simulation.

## Uniqueness
Depending on the intent of the invocation, you probably do not
want duplicate molecules generated. Uniqueness can be enforced
either within the molecules generated from a single starting
point, or across all molecules considered. Specify this with
either `-U each` or `-U all`.

## Output
Normal output just writes the changed molecules, but for exploration
and debugging, it can be useful to have the starting molecule written
at various times. The `-K` option controls output, and the following
qualifiers are recognised.

* **once**  parent is written at the first copy of each starting molecule
* **copy**  parent is written at the start of each copy.
* **every** parent is written before each variant

# Graph Edit Generation
Graph edit changes found in NextMove's SmallWorld are a subset of the
transformations found in `random_molecular_permutations`. Well, not
exactly the same, but somewhat close. Some of the transformations
built in here correspond to two edits in SmallWorld, for example
add an Oxygen atom, which I think would be two steps in SmallWorld -
add a Carbon, transform Carbon to Oxygen.

By creating
a probability file that assigns nonzero probability only to the
graph edit type transformations, we can mimic SmallWorld's
transformations - in a generative way. The tool is
```
LILLYMOL_HOME/contrib/bin/graph_edit_changes.sh
```
which uses `random_molecular_permutations` with a configuration
file that only enables graph edit changes. All options to the
tool are available.

Here are the nonzero probabilities from the config file.

```
add_double_bond 1.0
lower_bond_order 1.0
remove_atom 1.0
break_aliphatic_ring 1.0
make_ring 1.0
add_carbon 1.0
destroy_aromatic_ring 0.1
change_carbon_to_nitrogen 1.0
change_nitrogen_to_carbon 1.0
change_carbon_to_oxygen 1.0
change_oxygen_to_something 1.0
```

Some are very direct analogues of SmallWorld's graph edits,
add a carbon atom. But there is no remove bond directive,
since that would create multiple fragments - there are many
other directives that move things around, without creating
fragments.

## Diversion
This is not really documentation of the tool, just a record
of some exploration.

In exploring this here is an invocation that I found interesting
```
graph_edit_changes -U each -K every -x 1000 -p 1 -v file.smi
```
where for each molecules in the input, I make 1000 copies, and
from each of those 1000 copies, make a single graph edit change.

This took a couple of minutes, with the starting 2000 molecules
generating 84k new molecules, or about 42 variants per starting
molecule. The number of variants per molecule varied between 1
and 183. This took under two minutes.

The starting molecules look like
```
2000 molecules had between 7 and 50 atoms. Average 28.0095
```
The new molecules
```
84014 molecules had between 6 and 51 atoms. Average 29.9787
```
which looks not radically different - they have gained at
almost two atoms on average. How can that be, since these are
single edits? I think it is because the larger molecules
have more sites available, and so they bring up the numbers.

And because adding C, N, and O are separate transformations,
addition of atoms becomes more likely than a reduction.

We started with 2000 molecules, and asked for 1000 copies of
each, so we had 2M attempts at a transformation, but only
84k resulted in a unique variant.

If we look at how many atoms have a Hydrogen, 
```
new_tsubstructure.sh -s '[*H]' -v -v 
```
we find that on average, each starting molecule has about 8.6
atoms with an implicit Hydrogen. So if each of those sites were
replaced first with Carbon, then Nitrogen, then Oxygen, we
already have 26 variants per molecule before considering
other transformations.

So an average of 40 per starting molecule seems quite plausible.

But what this tells me is that the world of made molecules
is fairly sparse within the world of what could be made.
Now, many of these generated molecules could not be made, but
the point remains, occupied space is sparse.

# PostScript
New tool, `minor_changes` has been specifically developed to perform
simple changes. It is deterministic, unlike random_molecular_permutations,
and so can be more efficient in exploring structure space. This
distribution includes the minor_changes executable, but no fragment files.
In subsequent releases, fragments from public sources may be included.
