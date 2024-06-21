# JuLyMol

This is a work in progress.

Julia bindings work, but there is complexity in configuring CxxWrap
and other Julia dependencies.

As expected it is fast - usually much faster than python.

There are some unresolved issues around atom numbering and
iteration due to Julia's hyper-annoying choice of 1 indexing.
Atom numbers in LillyMol will always be zero indexed.
