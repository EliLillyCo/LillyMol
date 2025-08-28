# make_descriptors

This is a a descriptor generation tool that computes a variety of
molecular descriptors, using LillyMol tools.

> [!NOTE]
> Computation of 3D descriptors requires Corina in order to generate
> 3D structures from smiles. If you do not have corina, just omit
> all 3D descriptors. 

## Background
QSAR models can be built using LillyMol tools for descriptor computation.
`make_descriptors` is a convenient means of computing the various molecular
descriptors available in LillyMol.

The tool collects the set of descriptors to be computed from the command
line, and writes a Makefile that can evaluate those descriptors, and concatenate
the resulting files. Because make is used for coordinating the computations
parallelism is enabled using the parallelism built into make, which is very
efficient. The script supports a `-j` option, which is then passed to `make`.

INvoking make_descriptors with no arguments offers the following descriptor
sets available for computation
```
are recognised (* means 3D)
 -abr
 -ap
 -chgfp
 -cmi (vendor)
 -cnk     ** not yet working in LillyMol
 -comma*
 -dbf*
 -estate
 -ghose
 -ha 
 -hb
 -hpo
 -jurs*
 -jwc
 -jwdist*
 -jwdp*
 -jwmc
 -marvin (vendor)
 -medv
 -mk
 -morse*
 -pd     ** not yet working in LillyMol
 -sh*
 -tt
```
Again, note that descriptor sets marked with an asterisk are 3D descriptors which
depend on the tool (rcorina)[/src/Vendor/rcorina.cc] which depends on a shared
library from Corina, which must be separately licensed. If you do not have a Corina
license, just avoid all 3D descriptors.

If invoked with the -speed option, the script will report the relative
computational speeds of the different descriptor sets. If rcorina is available,
then all 3D descriptors will be by far the slowest to compute, since using Corina
to generate 3D structures is much slower than many of the other computations here.
Generally using Corina to generate 3D structures is very fast, but many of these
descriptor computations are faster.

A typical invocation might look like

```
make_descriptors.sh -w -abr -hpo -jwmc -j 4 train.smi > train.dat
```
which computes four different descriptor types, four way parallel.

If you wish to see the commands used, run as
```
make_descriptors.sh ... -o output.dat -v file.smi
```
and instead of stdout, 'output.dat' will be created and make will
echo the commands issued to stdout.

## Descrioptors
### abr
These are Abraham descriptors related to solubility and related effects
[abraham.cc](/src/Molecule_Tools/abraham.cc). These are often some of the most
beneficial descriptors to include. There are only 8 descriptors and they are
fast to compute.

### ap
Another fragment additive model from Michael Abraham. Fast to compute
and sometimes useful.

### chgfp
Apply the LillyMol formal charge assigner to a molecule and generate 2048 column linear
fingerprints that are sensitive to formal charge. 

### cmi
Fingerprint BioByte's clogp calculation. Reuquires a BioByte license.

### comma 3D

### dbf 3D
Distance between pharmacaphore features. Donors, acceptors and charged
atoms are defined and distances between pharmacaphore types are generated.
Also computes 2D values.
[estate](/src/Molecule_Tools/dbf.cc)

### estate
E-State descriptors.
[estate](/src/Molecule_Tools/jwestate.cc)

### ghose
The queries used in a Ghose Crippen logP estimation used as features.

### ha
Patterns of localised heteratoms.

### hb
Patterns of Hydrogen bonding features.

### hpo
Patterns of likely hydrophobic sections of a molecule.
[hpo](/src/Molecule_Tools/hydrophobic_sections.cc)

### jurs 3D
Surface area related features from Peter Jurs
[jwsa](/src/Molecule_Tools/jw_sa_db_descriptor.cc)

### jwc
CATS features [jwcats](/src/Molecule_Tools/jwcats.cc)

### jwdist 3D
Distribution of features along the longest spatial axis.
[maccskeys](/src/Molecule_Tools/jw_distribution_along_longest_axis.cc).

### jwdp 3D

### jwmc
Molecular connectivity indices.
[maccskeys](/src/Molecule_Tools/jw_molconn.cc).

### marvin (vendor)
Post process the output of cxcalc logP and logD.

### medv
MEDV descriptors [jwmedv](/src/Molecule_Tools/jw_MEDV.cc)

### mk
A variant on the original MACCS keys (192 features).
[maccskeys](/src/Molecule_Tools/maccskeys_fn5.cc).

### morse 3D

### pd
An expensive T-shaped fingerprint. Needs to be reimplemented.
[tshadow](/src/Molecule_Tools/iwpathd.cc).

### sh 3D
Shadow descriptors. Aligns and projects the molecule along axes.
[tshadow](/src/Molecule_Tools/tshadow.cc).

### tt
Topological Torsion derived features.

