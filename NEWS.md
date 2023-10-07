# June 2023

This release of LillyMol contains a great many updates.

Among the most notable are

# Build with bazel/bazelisk.
* Protocl Buffers used for substructure searching and reactions.
* Many additional substructure query functions.
* Many more utilities included in the release.

The build is now very different.  Previously LillyMol had tried hard
to ensure a lack of third party dependencies.  This became untenable
and now there are several external dependencies.  In addition,
building is now done with `bazel`.  This offers many advantages but
also makes building more complex than before.  See the README file for
instructions.

Switching to text format protocol buffers offers many advantages,
and will be used increasingly within LillyMol.
Old style query and reaction files will remain supported, but no new
functionality is being added to them. In order to gain access to new
functionality, you will need to use the text format protocol buffers.

Smarts provides a language for specifying matched atoms. But many
substructure query concepts are for groups of atoms, rings, substituents,
scaffolds, linkers... There are now more substructure abstract concepts
implemented, see especially the `global_id` concept which now
links the abstract concepts with matched atoms.

The no matched atoms between directive, `...` is now much more feature
complete, and allows specification of what atoms can, or cannot, be in
a region. Arbitrary element names are now more fully supported.

Within the code itself, new code is being implemented with a coding
style similar to what Google uses. clang-format is being used to 
transition parts of the code, and address sanitizer is used
during testing. GoogleTest is used for C++ unit tests - much more work
remains to be done there, but today there are over 25k lines of unit test.

## Future
As useful substructure search concepts are needed, they will be added
to LillyMol. If you have needs that are not addressed, please reach out.

A recent emphasis on LillyMol has been in the area of supporting and
augmenting AI driven generative designs. That will continue with
future work focussed on bringing precedented functional groups to
the generative space, combining the strengths of cheminformatics
enumeration and replacement based approaches and generative models.
Several tools supporting this are in this distribution.

