# ALOGP Models!

Overall I have found that this does not work well enough to be useful. Cool idea,
but not useful.
Do not use.


The LillyMol implementation of alogp was converted to read the parameters of the
model from a textproto. The tool has a default set of parameters from the paper,
and that default set is included in the file /data/alogp/default.textproto. A
parameter file like this can be provided to alogp via the -C option
```
alogp -C /path/to/config.textproto file.smi > file.alogp
```

alogp is relatively fast. It can score 59k molecules in just over 1 second.
This raises the possibility of using an optimisation tool to build a
fragment additivity model for any given objective - a docking score,
or other expensive computation, or an experimental measurement. The
optimisation algorithm would optimise the individual group contribution values
to best match the objective - the observed or computed values.
This has been done many times, there is nothing new here, just a convenient
way of doing it.

## ALOGP Models
The first task is to assemble the available data. If the data is experimental
you don't have a lot of choices to make - you still need a test set. But if the data is computational
we would suggest using [parsimonious_sample](/docs/Molecule_Tools/parsimonious_sample.md)
to identify a parsimonious sample of molecules that will constitute the
training set. Perform the expensive computation on that subset, as well as
a test set.

Even if the data is experimental, using parsimonious_sample might still be
a good idea.

## Building a Model
The Optimization.jl package from Julia is used to perform the optimization.
Optimisation is slow. Even though alogp itself can score a 50k model in about
one second, the model has 74+ adjustable parameters that must be
fitted. The Optimisation.jl package imposes a default number of iterations
of 10k, but this is likely much too few. Suggest 50k or more, depending
on the size of the dataset, and be prepared to wait overnight. At
one second per iteration, 50k seconds is about 14 hours.

Once the optimisation is complete, you can of course re-start
optimising using the final state produced by an earlier optimisation.
There will be no run-time penalty on the resulting model regardless
of how much optimisation is done.

Optimisation.jl points out that graident assisted optimisation can
be significantly beneficial. In the current implementation we are
using a no-gradient optimisation. Initial attempts to use an
automatic differentation option with Optimisation.jl indicated
significant extra complexity, so that was not pursued. It is also
unclear that in this case, a beneficial derivative could be computed.
Maybe, not sure.

For a new model, you will have no information about what parameters
might be useful, so I would suggest just use the optimised parameters
for alogp - your expermintal values might be related to logP anyway. 

### HOWTO
Three files are needed, all should be space separated.

1. A file of experimental values
2. A smiles file - same identifiers as are in the experimental value file.
3. An alogp::AlogpParameters textproto file containing the starting point.

Take a look at `${LILLYMOL_HOME}/data/alogp/default.textproto` for an example. These
are the default parameters for the alogP modelfor an example. As suggested
these are likely to be a reasonable starting point for any model
optimisation.
