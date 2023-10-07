# Smiles

Smiles in LillyMol are mostly standard, and can be consumed by many
other software packages.

## Aromaticity
Different software tools may make differing assessments of what is
aromatic. When presented with aromatic smiles from package `A`, package
`B` may be unable to parse the smiles, because package `B` do not believe
that this ring can be aromatic.

Also, fundamental to LillyMol is the need for a Kekule structure, so when
aromatic smiles are read, considerable effort must be undertaken in order
to find a valid Kekule form - which may occasionally fail.

We therefore prefer to store smiles as Kekule forms. That way we maximise
the probability that any other package can read the smiles. The resulting
extra size of the smiles is more than compensated for by the faster
processing enabled.

## Standardisation
There is a separate file on standardisation in LillyMol (link when done).
One of the decisions made in LillyMol is to prefer five valent Nitrogen atoms
over charge separated representations. Neither is correct, and they each
have advantages and disadvantages. Other packages may not be able to
consume five valent Nitrogen atoms. `fileconv -g rvnv5` will attempt
to convert all five valent Nitrogen atoms into charge separated forms
for use by third party software.

## Geometry
Coordinates can be added to a LillyMol smiles. For example benzene might be
```
C{{-0.0167,1.3781,0.0096}}1=C{{0.0021,-0.0041,0.002}}C{{1.2084,-0.6789,-0.0117}}=C{{2.3961,0.0285,-0.0201}}C{{2.3773,1.4107,-0.013}}=C{{1.1709,2.0855,0.0021}}1 benzene
```
where the atomic coordiantes are appended after each atom. This is usually
smaller than a .sdf file, and retains the convenience of having each structure
on one line. Tools that just read the connection table consume these smiles
just as if the coordinates were not there.

Over the years I have devoted considerable effort to try and find a compact
representation of the coodinates, with the objective of being able to write
the coordinates as something like `C{{121823}}=C{{121904}}` but such
efforts have generally not been successful, resulting in files that
are not that much smaller than the smiles file above. There are many
problems with this idea, and it has been abandoned, even though aspects
remain in the code. I have also tried working with binary proto files,
but again, the advantages are minor.

## Chemaxon Extensions
Chemaxon provides a wide variety of extensions to smiles via the `|...|` construct
after the smiles. LillyMol will recognise some of these, but most relate
to concepts that are not implemented in LillyMol - advanced stereochemistry, or
things that LillyMol can do via other means - arbitrary atomic symbols. Generally
support for Chemaxon extensions is slim... It will recognise coordinates.
