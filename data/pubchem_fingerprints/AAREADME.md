# Pubchem Fingerprints
This is an implementation of the Pubchem Substructure Fingerprint set
described at [ncbi](https://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.pdf).

While this implementation has not been extensively tested, it does appear to do a
pretty good job of implementing what is described in the description. But with so many
features described, errors are inevitable.

For tasks of model building and active retrieval, these fingerprints are generally
moderately useful. The script (pubchem_fingerprints.sh](/contrib/script/sh/pubchem_fingerprints.sh)
invokes the pubchem_fingerprints executable with these queries in place.
Computation speed is slow with 10k molecules being processed in 22 seconds - iwecfp
processes the same file in 1.5 seconds.
