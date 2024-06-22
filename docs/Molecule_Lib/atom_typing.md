# Atom Typing
Many tools within LillyMol use custom atom typing. With svmfp
models we have found that using more complex atom typing usually leads
to improved model performance. Over the years, the atom typing built
into LillyMol has become more complex and more fully featured.

It is important to note that as implemented, assigning atom types
and generating fingerprints are completely independent actions.
The fingerprint generator receives a set of atom types and has
no knowledge of what they might represent. So, whether these atom
types get used for an EC fingerprint, linear fingerprints, or
atom pair fingerprints is arbitrary. So, whether these atom
types get used for an EC fingerprint, linear fingerprints, or
atom pair fingerprints is arbitrary.

## History
Once the idea that atom types can be more complex entities than
just the atomic number, there is essentially 
no limit to what can be done with atom typing. This raises both
many opportunities and many challenges.

Initially we started building combinations of atomic properties and
giving those combinations names. Today, the following names are still
recognised

- bs
- cc
- ch
- complex
- expt
- none
- nox
- pp
- sb
- sf
- sfx
- tt

Some of these correspond to composites of individual atomic
properties, and could also be specified via the new `UST:` (User Specified
Type) directive described below.  Others are unique.  Clearly
assigning a unique text name to every individual combination of atom
attributes was not going to work.  Therefore the `UST:` directive was
added.

The named atom types are

### bs
So-called basic atom type. There are 6 different atom types.
1. aromatic atoms.
2. singly connected carbon atoms
3. all other carbon atoms
4. fully saturated, has hydrogen attached
5. fully saturated, no hydrogens
6. with a hydrogen
7. other

This crude atom type sometimes works well in models.

### cc
First applies the `complex` atom type. If the molecule has no formal charges
return that. For every atom that has a formal charge, perturb the atom type
depending on whether it is positive or negative.

Since chemical standardisation tries hard to drive out formal charges,
it is not clear that this will be useful. 

### ch
Carbon atoms get 1 type, non carbons get a different type.

Very simplistic, but sometimes useful in similarity searches.

### complex
One of our best atom types, especially in models. Ultimately it is
```
  invariant[i] = 100 * hcount + 10 * pi + 4 * arom + m.ncon(i);
```
which is a function if the number of hydrogens attached to an atom,
the pi electrons, whether or not the atom is aromatic, and the number
of connections. Since almost all calculations are done with implicit
rather than explicit Hydrogens, `hcount` will be implicit Hydrogens.

Note that the element type is not part of the atom typing, although it
is quite possible that due to the constraints imposed by valences, it
might be implicitly encoded.

### expt
An *experimental* atom type. This is for developers to try different
things, which if successful might then get their own name.

### none
All atoms get the same type. This is a step towards fingerprinting
the molecule as a skeleton. But note that this atom type just controls
the atom type, it does not say anything about bonding.

### nox
Consists of 3 atom types

1. carbon
2. oxygen or nitrogen
3. other

This idea proved useful in a specific project once upon a time.

### pp
This is the pharmacaphore fingerprint. This depends on externally
specified configurations for formal charge and Hydrogen Bonds. It is
complex, and is described in its own section below.

### sb
This is an atom typing designed to mimic the atom typing found in
Sybyl, an extinct molecular modeling package. It is seldom useful.

### sf
This is a type designed to capture atomic properties relevant to
synthetic precedent databases. After some experimentation, the `sfx`
typing was found to work better.

This type consists of

1. atomic number
2. aromaticity
3. ring or not
4. saturation

### sfx
The atom typing currently used for synthetic precedent databases.
It consists of the following attributes.

1. ring status - there is a complex algorithm that looks at ring fusion, 
aromatic or aliphatic.
2. atomic number
3. number of connections
4. hydrogen count
5. formal charge

### tt
An atom typing derived from the original Topological Torsion paper. It
consists of

1. atomic number
2. number of connections
3. hydrogen count

This very simple atom typing often does quite well in models.

### Pharmacaphore Atom Type
This is driven by a configuration file that contains directives on how
to form a pharmacaphore fingerprint. Those directives wil include

- charge assigner
- donor acceptor assigner
- hydrophobe queries
- atom type for non pharmacaphore atoms
- combine donor and acceptor

The last directive, if present, says do not differentiate donors and
acceptors, they all get the same type. 

Note that there are some fundamental problems with this atom typing. Some
atoms can be both donors and acceptors, but atom typing depends on 
there being a unique type for an atom. By default, a donor will be
a different atom type from a dual, which seems undesirable. The only
way around this would probably be to fingerprint the molecule multiple
times, once with all duals as donors and once with them all as
acceptors. But fingerprint generators do not know anything about this
so that is not really feasible.

Importantly, specifying a type for the non pharmacaphore atoms is
very important. This derives from work at Glaxo I recall where their
efforts at implementing pharmacaphore atom types was floundering. And
then they introduced an 'other' type, and things started working.
This is hardly surprising since pharmacaphoric features usually
comprise only a small fraction of the atoms in a molecule.

This atom typing seldom works well, probably because of the limitations
noted above. In order to be useful, a custom fingerprint generator
would likely be needed - something we have never implemented.

## UST
The User Specified Type directive allows combination of any number of
atomic properties - subject to some constraints. As an implementation
detail, arbitrary numbers are assigned for each atom type, and
these are multiplied and added - basically in the hope that there will not be
too many collisions. Generally this works OK.

The following atomic properties are recognised with a UST atom typing

- A aromatic or not.
- C number of connections.
- E heteroatom or not.
- F for ring atoms, is the ring isolated or fused.
- G element type - useful for atoms that do not have an atomic number. *
- H number of hydrogens attached.
- I isotope.
- K connected atoms hash. Combines the atomic numbers of neighbours.
- L size of largest ring.
- M first assign `Y` types, then all aromatic atoms become identical. *
- N none - all atoms get the same type. *
- O formal charge.
- P pi electrons - counts differentiated.
- Q pi electrons - present or absent.
- R number of ring bonds attached.
- S size of smallest ring.
- T tautomer inspired atom type. likely tautomeric nitrogens all get the same type. *
- U unsaturation - includes aromatic.
- X atom type adjusted for `centrality` of atom in the molecule.
- Y compressed atomic number *
- Z atomic number *

Atomic properties marked with an asterisk reset the atom type, whereas
the non asterisk properties alter an existing value. In the code, those
that reset the atom type are run first. But note that combining properties
that each reset the value will result in unpredictable results, without
warning. `UST:YZ` is undefined, although it will silently yield a valid
result.

Some of these are clearly better ideas than others. Properties `F`, `K` and
`X` will likely be repurposed, so do not use.

Common atom types might include things like

- UST:ACHPY
- UST:ACRY
- UST:ACHPRT

by convention these are alphabetized, which will hopefully lessen the
probability of unfortunate letter combinations.

The `Y` atomic property is interesting. I have never seen a qsar model
that was improved by differentiating Cl, Br and I atoms. So, the `Y`
atom type compresses all those types to Cl. Always prefer this over `Z`
unless you have a specific reason.

The `T` atomic property assigns the same atom types to aromatic nitrogen
atoms in the same ring, where one might need to be drawn with a Hydrogen
and the other without, while in reality they might be tautomerically
equivalent. Certain kinds of Oxygen and Sulphur atoms are made
equivalent, and from an aromatic ring, singly connected oxygen
atoms are merged. This atom typing often works well.

Note that when used in gfp_make, the `UST:...` fingerprint types are specified
as things like
```
-EC3:AHY
-EC3:ART
-MAP8:CPY1
-IW7:AY
```
which gfp_make will expand to
```
... -P UST:...
```
during the invocation.

## Shell Expansion
All `UST:` directives can be followed with a digit. When this is
specified, atom types are assigned. Then atom types are re-assigned
based on the atom types of the atom and its neighbours - just like
what happens when EC fingerprints are formed. Obviously this would
not make sense if EC type fingerprints were being generated, but
we have seen good results when using atom pairs - in effect, we are
fingerprinting dumbell type shapes. I have only ever used 1 level of shell
expansion.

## Implementation

Within the `c++` code, there is an Atom_Typing object that usually
does two things within a program.

1. Gets initialised via the command line
2. For each molecule read, assign the atom types.

In pseudocode that often looks like

```
// If the -P option is given, pass that to build.
if (cl.option_present('P'))
  if (! atom_typing.build(cl.value('P'))) {
...
// If it has been activated, assign atom types.
std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);
if atom_typing.active()
  atom_typing.process(molecule, atype.get())
```

Because this object can be compiled into many different tools,
we get the same behaviour across many different tools. Generally
we have tried to make the `-P` option be used for atom
typing specifications, so a typical invocation might look like
```
tool <other options> -P UST:ACY1 input.smi > output
```

## External Query Driven
For maximum flexibility the atom typing allows specification of
functional groups via an `atom_typing_spec::External` proto.

That proto is very simple, consisting of repeated instances of

```
message QueryAndValue {
  // A query type that defines one or more queries.
  // SMARTS:C for example.
  repeated string query = 1;

  // The atom type value assigned to atoms that match
  // the query.
  optional uint32 value = 2;

  // Ignored
  repeated string comment = 3;
}
```

For each molecular feature, there is a query specification, and a
numeric value that is assigned to atoms covered by the query. See
the [file](external.textproto) for an example.

This can be applied as
```
-P EXT:/path/to/external.textproto
```

## Future
As is often the case, hindsight is very good and there were design
decisions made many years ago that today are not so great. Work
continues on finding useful atomic descriptions and their combinations
that demonstrate utility in models.

Note that it is usually unpredicable which atom type combinations
will work best in any given situation.
