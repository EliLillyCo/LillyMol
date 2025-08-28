# Fingerprints
Long term, the LillyMol GFP (Generalised FingerPrints) class will be
made available from python. In the meantime there are some capabilities
implemented today that may be useful for computing small numbers of similarities.

## Similarity
There is no best similarity measure. The best similarity measure would
likely be one that accurately mapped to changes in Biological activity.
Much of Cheminformatics, and other, research is designed to find
means of accurately predicting similarity in Biological activity.
It is hard.

A good fingerprint might be the one that does the best at tracking
changes in Biological activity. That will likely be target dependent
and will need extensive study to identify.

A good fingerprint might be one that generally corresponds with human
perceptions of chemical similarity.

Our experience is that fingerprints like EC (Extended Connectivity, or Morgan)
fingerprints work best for things like SVM fingerprint models, linear
path fingerprints tend to work best for corresponding to human
perception. It all depends.

## HowTo

The following toy application shows a simple N*N near neighbour
computation in LillyMol python.
```
from absl import app
from absl import logging

from lillymol import *
from lillymol_io import *
from lillymol_fingerprint import *

def main(argv):
  if len(argv) == 1:
    logging.error("Must specify input file")

  mols = slurp(argv[1])
  logging.info("Read %d molecules from %s", len(mols), argv[1])

  fps = [linear_fingerprint(mol) for mol in mols]
  logging.info("Fingerprints generated")

  nfp = len(fps);
  for i in range(nfp):
    max_similarity = 0.0
    idmax = -1
    for j in range(nfp):
      if i == j:
        continue
      t = tanimoto(fps[i], fps[j])
      if t > max_similarity:
        max_similarity = t
        idmax = j

    print(f"{mols[i].smiles()} {mols[i].name()} {mols[idmax].smiles()} {mols[idmax].name()} {1.0 - max_similarity}")

if __name__ == "__main__":
  app.run(main)
```
This does a very dumb N*N nearest neigbour computation. With some book-keeping
only half the calculations need to be performed. This example is just designed
to demonstrate outlines of how fingerprints in LillyMol python work.

Unfortunately this is quite slow. Running 2000 random Chembl molecules takes 21 seconds.
Running gfp_nearneighbours_single_file on the same input file, but with standard
gfp fingerprints, takes less than 1 second. But that is dealing with bit vector
fingerprints and can use popc instructions for computing similarity. This calculation
is using counted fingerprints, which will necessarily be a much more expensive
computation.

When GFP fingerprints are ported into the LillyMol Python environment, faster
computation of binary fingerprints will become available. That said, there
are advantages to counted fingerprints, since they do not suffer from the
repeated feature problem of binary fingerprints. If a molecule contains
four instances of a feature, that will still just set the bit once in a binary
fingerprint, whereas in a counted fingerprint, the number of instances will
be recorded and will count in the similarity computation.

Only linear fingerprints have a method for constructing fingerprints such as the above.
All fingerprints can work via a fingerprint generator, that can be configured with
things like atom typing.
```
  fpgen = LinearFingerprintCreator(2048)
  fpgen.set_max_length(7)
  fpgen.set_atom_type("UST:AY")

  fps = [fpgen.fingerprint(mol) for mol in mols]
  logging.info("Fingerprints generated")
```

For example to generate 1024 bit EC fingerprints, of radius 3, using an atom type of 'UST:AY' that
could be done by
```
  fpgen = ECFingerprintCreator(1024)
  fpgen.set_max_radius(3)
  fpgen.set_atom_type("UST:AY")

  fps = [fpgen.fingerprint(mol) for mol in mols]
```

Atom pair fingerprints are similar, this time restricting separations to 10 bonds.
```
  fpgen = AtomPairFingerprintCreator(1024)
  fpget.set_max_separation(10)
  fpgen.set_atom_type("UST:AY1")

  fps = [fpgen.fingerprint(mol) for mol in mols]
```
All generate numpy byte arrays, containing counted fingerprints.

Again, speed is not good. This can be improved by several approaches

* bitvector based fingerprints
* GFP fingerprints
* A multiple fingerprint container that can be queried with just 1 call from python.

And all constructors could be altered to allow pythonic **kwargs to deal with
settable parameters.
