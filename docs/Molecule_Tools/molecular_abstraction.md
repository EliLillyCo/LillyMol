# Molecular Abstractions
Puts a formal charge of 1 on every atom.
`molecular_abstractions` is a tool which converts an input molecule
to various changed and abstract forms.

## HOWTO
The changes are specified via the `-a` option, and consist of
directives, which may contain extra instructions. The embedded
instructions might specify information about what to change, or
they might direct that the molecule, in its current form, be
written or fingerprinted.

An example usage might be
```
molecular_abstractions -a 'scaffold(WRITE)' file.smi
```
which converts the input molecule to the scaffold form, and
since the `scaffold` directive has a `WRITE` directive, a
the smiles of the scaffold is produced.

To generate the fingerprint of the scaffold,
```
molecular_abstractions -a 'scaffold(FP)' file.smi
```

Directives can be chained, so to generate a molecular
skepeton, where all atoms are carbon, all bonds
are single bonds, and all formal charge values set
to zero.
```
molecular_abstraction.sh  -a 'allatoms.allbonds.charge(WRITE)' file.smi
```

The empty `allatoms` directive is a shortcut for setting all atoms
to a single type, which is `C` by default. Same with `allbonds`
where the default behaviour is to convert all bonds to single bonds.
The default behaviour for the `charge` directive is to place
zero formal charge on all atoms.

Clearly all these particular directives could be made more
efficient by specifying changes only to those atoms that needed
to be changed. See below for how to more finely tune these
actions.

## Interesting cases
Many operators also recognise the `ISO` directive, which means apply
isotopic labels to the resulting molecule.  Where the isotopes are
placed, and whether isotopes make sense for an operator, depends on
the operation.

So, if you wanted the largest ring system, labelled by attachment point, that would be:
```
-a 'bigring(ISO WRITE)'
```
This applies a default isotope of 1 to the attachment points.  Use
`ISO=nn` to specify your own isotopic label.  If you then wanted all
aliphatic heteroatoms transformed to Sulphur
```
-a 'bigring(ISO).translate([!#6;R]=S WRITE)'
```
If you wanted all rings, but without doubly bonded =O atoms:
```
-a 'rings.rmat([D1]=[R] WRITE)'
```
That should do it. Note that only the first matched atom is removed.

By default, atoms are simply excised, with no rejoining.  The `rejoin`
directive tries to rejoin the connections lost.  So, to convert esters
to ketones try:
```
-a 'rmat([OD2R0]-C=O rejoin WRITE)'
```

Atoms that have higher connectivity can be rejoined with the `RJ=n` directive.
So, if you'd like to rejoin items that had 3 connections, try:
```
-a 'rmat(smarts RJ=3 WRITE)'
```
To remove multiple atoms, just specify the smarts:
```
-a 'rmatoms(smarts)'
```
So, to remove ring attached cyclopropyl rings:
```
-a 'rmatoms([/IWxR]-!@C1-[CD2]-[CD2]1)'
```
The construct `[/IWxR]` is interpreted as an atom in a Ring.  The `/IWx`
directive means exclude this atom from the matches reported, so the
ring atom will not be removed. Recursive smarts can also be used
to exclude atoms from what is reported as matched.


To specify a changed bond type, you need to specify a smarts that
defines two matched atoms.  To convert all aliphatic, ring double
bonds to single bonds, try
```
-a 'cbt([A;R]@=[A;R]=- WRITE)'
```
To remove a bond, specify ‘.’ as the replacement bond. So, to remove
all single bonds between rings: 
```
-a 'cbt([R]-!@[R]=. WRITE)'
```
You can change all bonds with the allbonds directive, and all atoms with the allatoms directive. To 
transform to the skeleton, try:
```
-a 'allbonds.allatoms(C WRITE)'
```
## Directives
The directives currently recognised are:
* **allatoms** Change all atom types 
* **allbonds** Change all bonds to a given type (default single)
* **arf** Abstract Ring Form 
* **bigring** Remove all but the largest ring system
* **cbt** Change bond types between matched atoms
* **charge** place formal charges 
* **comprconsec** compress consecutive CH2 groups to just one 
* **frag** Fragment removal and selection 
* **invscaf** The atoms that are **not** part of the scaffold.
* **isotope** place isotopes 
* **rings** Remove all non-ring atoms 
* **rmat** Remove single atom 
* **rmatoms** Remove multiple atoms 
* **rmbond** remove bond 
* **rmrd2** remove all two connected ring atoms - compresses rings
* **rplink** Replace inter-ring Linker atoms 
* **scaffold** Form the molecular scaffold
* **sss** filter by substructure 
* **translate** Change specified atom types 

### allatoms
Converts all atoms to a single type. Change all atoms to Sulphur
```
-a 'allatoms(S WRITE)'
```
bonding and charges not altered. If you only want to change certain atoms
use the `translate` directive.

If no element is specified, all atoms are converted to Carbon.

### allbonds
Change all bonds in the molecule to a single type.
```
-a 'allbonds(WRITE)'
```
The default is that all bonds become single bonds. To change all bonds
to double bonds
```
-a 'allbonds(= WRITE)'
```
which is unlikely to be useful. If you want to change specific bonds,
use the `cbt` directive.

### arf
Abstract Ring Forms.This transformation replaces aromatic rings
with a single atom, `Ar`, aliphatic rings with `Al`. If the rings
are fused, they will be connected by an `=` bond, and if strongly
fused, alphatic only, via a triple bond.
You can change the elements applied with the 
`AROM=El` and `ALIPH=El` directives. Spiro fusions are split apart
and a dummy `Sg` atom inserted.

The following qualifiers are recognised.
* **lrs** label the abstract atoms with the ring size.
* **lhc** label the abstract atoms with the heteromat count in the original ring.
```
-a 'arf(lrs WRITE)'
```
Note that also specifying an ISO= directive would be undefined - should be
detected as an error. TODO: ianwatson

For example with input `C1CC1Cc1ccccc1` running
```
molecular_abstraction.sh -a 'arf(WRITE)'
generates
[Ar]C[Al]  ARF
```
Benzimidazole `c12ccccc1ncn2` produces `[Ar]=[Ar]` and adamantane,
`C1C2CC3CC1CC(C2)C3` generates `[Al]1#[Al]#[Al]#1` where a ring
structure is also used. Using `-a 'arf(lrs WRITE)'` on Benzimidazole
results in `[6Ar]=[5Ar]` where the isotope indicates the size of the ring.

These reduced graph forms can be further optimized with other
transformations to build custom abstract forms.

## bigring
Removes all atoms except those in the largest ring system. Add the spiro
directive to have ring systems span across spiro fusions. Specify an
isotope and the join points will be labelled.
```
-a 'bigring(spiro ISO=3 WRITE)'
```

## cbt
Change the bond type between two matched atoms. In the smarts, the two
atoms defining the bond are assumed to be the first two matched atoms.
```
-a 'cbt(S=O->- WRITE)'
```
Note that `=-` also works, but seems awkward.

## charge
Set or alter a formal charge value.

If no directives are present, all charges are removed.
```
-a 'charge(WRITE)'
```

If just one token is present, it
is assumed to be a smarts, and a charge of 0 is placed on the matched
atoms. The neutralize acids
```
-a 'charge([O-]-C=O WRITE)'
```
To set the a charge on a specific matched atom
```
-a 'charge(O=C-[OH] 2=-1 WRITE)'
```
would place a negative charge on matched atom 2, the OH. That is equivalent
to
```
-a 'charge(SMARTS=O=C-[OH] 2=-1 WRITE)'
```
which might be clearer.
```
-a 'charge(1 WRITE)'
```
Puts a formal charge of 1 on every atom.

### comprconsec
Change all consecutive -CH2- groups to just one CH2.
```
-a 'comprconsec(WRITE)'
```
If you add the
ISO directive, the remaining CH2 group will get an isotopic label
corresponding to the initial length of the chain.
```
-a 'comprconsec(ISO WRITE)'
```

If you want some other atom type, just use a smarts.
```
-a 'comprconsec([CD2R1H2] WRITE)'
```
This will compress consecutive ring CH2 groups.

If you want to remove and compress all ring CH2 groups, use the `rmrd2`
directive which is specially formulated for the task.
Whereas `comprconsec` will always leave at least one CH2 
group, rmrd2 may remove all of them - leaving the substituted ring skeleton.

### frag
This is for controlling selection of fragments. Note that `fileconv`
has extensive capabilities for handling fragments, and that is
probably the most flexible fragment selection tool. Clearly it would
be desirable to combine the two and have an integrated fragment
handling functionality.

By default, all but the largest fragment are removed.
```
-a 'frag(WRITE)'
```
Specific fragments can be removed or kept.
To remove all fragments that match SMARTS
```
-a 'frag(REMOVE=SMARTS WRITE)'
```
and to keep only those fragments that match SMARTS
```
-a 'frag(KEEP=SMARTS WRITE)' 
```

### invscaf
Currently this appears broken. It is writing an interesting output,
but it is **not** the inverse scaffold. TODO ianwatson: fix.

The atoms that are **not** part of the saffold, also called
molecular spinach - the stuff hanging off the core.

```
-a 'invscaf(WRITE)'
```

### isotope
Sets isotopic labels. Note that `fileconv`
has extensive capabilities for handling isotopes, and that is
probably the most flexible isotope control tool. Clearly it would
be desirable to combine the two and have an integrated isotope
handling functionality.

Note that there is ambiguity between the `ISO=nn` that is recognised
by many directives, and this one. For clarity, do **not** use the
`ISO=` qualifier here - although there is an attempt to make it work.

The default behaviour is to remove all isotopes.
```
-a 'isotope(WRITE)'
```

Somewhat counterintuitively, the default if a query is specified is
to place isotope 1 on the matched atoms.
```
-a 'isotope(C WRITE)'
```
To place a different isotope, and on a different atom
```
-a 'isotope(SMARTS=[NH]-c 1=3 WRITE)'
```
which places isotope 3 on matched atom 1, where matched atom 1
is part of an aniline.

### rings
remove all non-ring atoms. Note that some doubly bonded exocyclic
connections that are needed to preserve aromaticity will also be
retained. There is an attempt to **not** destroy aromaticity.
```
-a 'rings(WRITE)'
```

### rmat
Remove a single atom. At a minimum a smarts must be specified.
```
-a 'rmat(c WRITE)'
```
will remove all aromatic carbon atoms!

The `rejoin` qualifier can be used to make the removal an
excision, and try to preserve the existing bonding. So if a
two connected atom is removed, the adjacent bonds will form
a bond. But note that if adjacent atoms are removed, it is
not smart enough to figure out that a bond could still be
made. This could be fixed.

For example if the molecule `NCCCO` is processed as
```
-a 'rmat([CD2] WRITE)'
```
the result is `O.N` which is correct, all 2 connected Carbon atoms
are removed. Note that it did not figure out that a join was possible.
But if instead we only remove the central carbon
```
-a 'rmat([CD2T0] rejoin WRITE)'
```
then the result is `OCCN`.

Note that there is an interesting feature if a three connected
atom is removed, and a rejoin is requested. It forms a three
membered ring.
So if `CC(C)C` is processed with
```
-a 'rmat([CD3] rejoin WRITE)'
```
the result is `C1CC1`. Not sure if this is a bug or a feature.

I am not sure why this directive was necessary.
It may be necessary to only attempt rejoining when the matched atom
has a certain number of connections. Specify one or more connectivity
values via the `RJ=n` directive.

### rmatoms
Used for removing many atoms. By default, all atoms matched by
the query are removed.
To remove carboxyllic acids
```
-a 'rmatoms([OH]-C=O WRITE)'
```

### rmbond
Remove the bond between two matched query atoms.
To break all biphenyl bonds
```
-a 'rmbond(a-!@a WRITE)'
```
You can also discard the fragment containing matched atom 0 or 
matched atom 1 with the REMOVE=0 or REMOVE=1 directives.
```
-a “rmbond([CD2T0]-[CD3]=O REMOVE=0 WRITE)”
```
This will break esters bonds and remove the OH side, leaving an aldehyde.

### rmrd2
This removes two connected ring atoms.
Convert `CC1C(C)CCC1` into `CC1CC1C` with
```
-a 'rmrd2(WRITE)'
```
To restrict to just Carbon atoms

```
-a 'rmrd2(C WRITE)'
```
Note that the extra token is interpreted as an atomic symbol, and
not as a smarts. That should probably be changed.

This directive is useful forming abstract ring forms.

## rplink
Replace linker atoms.

Linker atoms are atoms that are between rings. So, `C1CC1CCC1CC1` becomes
`[Li](C1CC1)C1CC1`. `C1CC1CCCCCCCCC1CC1` is also transformed to the same.
`C1CC1CCCC(C1CC1)CCCCC1CC1` is transformed into `C1C([Li](C2CC2)C2CC2)C1`.
```
-a 'rplink(WRITE)'
```

If an element other than Li is desired, just add another token
```
-a 'rplink(U WRITE)'
```

### Scaffold
The scaffold atoms. Ring atoms, and atoms between the rings, as well
as some doubly bonded atoms.
```
-a 'scaffold(WRITE)'
```

The `keepfirst` qualifier retains atoms that are directly attached to
scaffold atoms.

### sss
Substructure search. Molecules not matching the query have all their
atoms removed, and when that happens, processing stops.
```
-a 'sss(SMARTS=O=C-N WRITE)'
```
results in only molecules that contain an amide. Any number of `SMARTS=`
directives can be specified, and by default they function as an `or` match.
```
-a 'sss(SMARTS=O=C-N SMARTS=F WRITE)'
```
which is molecules that contain either an amide or a Fluorine.

The `NONM` qualifier inverts the meaning of the matching, so non
matching molecules are propagated.

Tokens not otherwise recognised are interpreted as the path to a
query file - anything recognised by the `-q` option to many tools.
So, a file of smarts could be specified as `S:/path/to/file`.

### translate
Translate individual atoms from one element to another.

Takes one or more directives of the forms `smarts=ele` or `smarts->ele`.
To translate all amides into esters
```
-a 'translate([ND1,ND2]-C=O->O WRITE)'
```
Note that it is the first query atom that is transformed. Note that
substructure perception is done after each change, so the quite
redundant transformation
```
-a 'translate(F->[U] [U]->Cl WRITE)'
```
works the same as `F->Cl`.

## Options.
The usage message is
```
Generates molecular abstractions
  -a ...        specify abstraction(s) to be created, enter '-a help for info
  -B <fname>    specify abstraction(s) in a file, same syntax
  -p            write the parent molecule
  -z ...        options for what to do when no changes during a stage
  -C            remove invalid chirality from output molecules
  -h            unfix any explicit implicit hydrogen specifications
  -c            remove all chirality from input molecules
  -l            reduce to largest fragment
  -t            remove cis-trans bonds from input
  -Y            standard fingerprint options, enter '-Y help' for info
  -f            work as a TDT filter
  -F <stem>     write final molecule to file <stem>
  -o <type>     type for -F file
  -i <type>     input specification
  -g ...        chemical standardisation options
  -E ...        standard element specifications
  -A ...        standard aromaticity specifications
  -K ...        standard smiles options
  -v            verbose output

```

The `-a` option has been explained above.

The `-B` option allows placing the -a directive in a file, avoiding
lengthy command lines and difficulties with quotes.

The `-p` option enables writing the parent molecule as well as
the transformed form.

The `-F` option allows specifying a destination other than stdout.

The `-z` option can control what happens when a change happens or
does not happen at a given stage. Enter `-z help` for more details.

The `-Y` option allows specifying different behaviour for any fingerprint
opterations. Enter `-Y help` for info.

During transformations, unusual valences can be encountered. The `-h`
option can sometime alleviate these by converting atoms with a known
hydrogen count, to a free hydrogen count.

The `-c` option removes chirality from input molecules.

During processing, existing chirality may be invalidated. Use the `-C`
option to re-examine chirality in the transformed molecules, and remove
any that is no longer valid.

The `-t` option removes any cis-trans bonding in the input. Recommended.

The `-l` option strips input molecules to the largest fragment.


