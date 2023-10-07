# LillyMol
## Value Proposition
LillyMol is a Cheminformatics toolset designed to support
drug discovery. It contains tools that have been found to be
useful in drug discovery, together with an underlying C++ library
all made available via a permissive license.

Being C++ LillyMol is usually fast.

LillyMol does not attempt to address all the important tasks
needed for Cheminformatics. It was developed in a heterogenous
environment where a great variety of open source and commercial
tools were available. Where good or excellent tools already
existed, LillyMol made no attempt to replicate that.

LillyMol is a command line, file based processing set. It
works very naturally with system parallel processing tools
since many large scale Cheminformatics problems are quite
`pleasingly parallel` and are amenable to processing in
sharded form. It also works well with pipelined command
processing, where the output from one command is fed to
the next command, thereby taking advantage of the multi-core
environments that today are uniquitous.

Since all LillyMol tools are built on a common code base,
tools tend to have common arguments and behaviours.

Integration with languages such as Julia and Python is possible, and
prototypical implementations have been done. As time allows
work will continue on that area.

The main value drivers for LillyMol include

* Many tools designed for drug discovery.
* New concepts for substructure searching.
* Novel means of reaction enumeration.
* A very flexible environment for similarity searching.
* High performance similarity searching.

## I/O
LillyMol uses common I/O functionality, see ... for details.
